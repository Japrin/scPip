#!/bin/bash

g_prj_id="BRCA"

cDir=`pwd`
### columns required in $file_obj_in: aid measurement platform seufile scefile resolution
### columns required in seurat or sce object: "patient"
file_obj_in="$cDir/list/obj.T.wRef.list"
file_int_in="$cDir/list/obj.T.wRef.list.r1.list"
file_limma_in="$cDir/list/obj.T.list.forLimma.list"
file_limma_out="$cDir/list/obj.T.list.limma.sc.list"
cellType="T"
oDir="$cDir/OUT.byDataset"
shDir="./sh.byDataset.$cellType"
mkdir -p $shDir

sDir=`R --slave -e 'sDir <- system.file("script",package="scPip"); cat(sDir)'`
echo $sDir

######### 1.1 generate scripts to run Seurat for each dataset
while read aid measurement platform seufile scefile resolution
do
(
cat <<-HERE
#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.$cellType.$aid.%j.out
#SBATCH -e S.$cellType.$aid.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/run.seurat3.basic.R \\
	-a $seufile \\
	-b $scefile \\
	-o $oDir/$cellType/$cellType.$aid/$cellType.$aid \\
	-d 15 \\
	-n 12 \\
    --keep globalC:T,T/NK \\
    --removeContamination plasmaB:0.75,caf:0.75,epi:0.75,mac:0.75 \\
    --deg \\
	--resolution $resolution \\
	-m $measurement \\
	--platform $platform
HERE
)>$shDir/byDataset.$cellType.$aid.sh
done< <(awk 'NR>1' $file_obj_in)

########## 1.2 submit jobs to run the scripts ########
#cd $shDir
#ls *.sh | awk '{print "sbatch "$0}' | bash
#cd $cDir
##################################################

######## 2.1 prepare file list for integration ############
sed '1,1d' $file_obj_in \
    | perl -ane 'BEGIN{ print "data.id\tmeasurement\tplatform\tdefile\tscefile\tseufile\n"  }
                chomp; print join("\t",@F[0..2],
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/limma/'$cellType'.$F[0].de.out.limma.rda",
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/'$cellType'.$F[0].sce.rds",
                "'$oDir'/'$cellType'/'$cellType'.$F[0]/'$cellType'.$F[0].seu.rds")."\n" ' \
    > $file_int_in
########################################################


######################## 2.2 first run ########################
(
cat <<-HERE
#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o S.int.$cellType.%j.out
#SBATCH -e S.int.$cellType.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/wrapper.run.inte.R \\
    --inFile $file_int_in \\
    --outPrefix $cDir/OUT.int.$cellType/int.$cellType \\
    --corVar S.Score,G2M.Score,DIG.Score1
HERE
)>$shDir/inte.$cellType.sh
########################

########## 2.3 submit jobs to run the scripts ########
#cd $shDir
#sbatch inte.$cellType.sh
#cd $cDir
##################################################

######################## 3.1 examine the result, find whether there are cells needed to be excluded
./w.checkContamination.T.R
########################

######################## 3.2 second run ########################
### if need to filter out some cells, add:
###    --excludeCells "$cDir/OUT.int.$cellType/int.$cellType.contamination.vec.rds" \\
(
cat <<-HERE
#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -o S.int.$cellType.2nd.%j.out
#SBATCH -e S.int.$cellType.2nd.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
$sDir/wrapper.run.inte.R \\
    --inFile $file_int_in \\
    --outPrefix $cDir/OUT.int.$cellType.2nd/int.$cellType \\
    --excludeCells "$cDir/OUT.int.$cellType/int.$cellType.contamination.vec.rds" \\
    --corVar S.Score,G2M.Score,DIG.Score1,ISG.Score1
HERE
)>$shDir/inte.$cellType.2nd.sh

########## 3.3 submit jobs to run the scripts ########
#cd $shDir
#sbatch inte.$cellType.2nd.sh
#cd $cDir
##################################################


######################## 4. cluster annotation ########################
./w.ann.$cellType.R

######################## 5.1 run limma per dataset ########################
join -1 1 -2 1 \
    <(cut -f 1-3 $file_int_in|sort -k 1r,1) \
    <(ls $cDir/OUT.int.$cellType.2nd/sce/*.sce.rds | perl -ane 'BEGIN{print "data.id\tscefile\n" } chomp; /sce\/(.+?).sce.rds/;print "$1\t$_\n"' | sort -k 1r,1) \
    | sed 's/\s\+/\t/g' \
    | awk -F"\t" -v OFS="\t" 'BEGIN{ print("data.id\tmeasurement\tplatform\tscefile") } !/^data.id/{print $0}' \
    > $file_limma_in

while read data_id measurement platform scefile
do
(
cat <<-HERE
#!/bin/bash
#SBATCH -p all
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -o S.limma.$cellType.$data_id.%j.out
#SBATCH -e S.limma.$cellType.$data_id.%j.err
#SBATCH --no-requeue
echo begin at: \`date\` host: \`hostname\`
echo \`hostname\`
$sDir/wrapper.run.limma.R \\
        -b $scefile \\
        -o $cDir/OUT.int.$cellType.2nd/limma.sc/$data_id/limma.sc.$data_id \\
        --platform $platform \\
        --group "majorCluster" \\
        --groupMode "multiAsTwo" \\
        -n 8 \\
        -m $measurement
HERE
)>$shDir/limma.$cellType.$data_id.sh
done < <(awk '!/^data.id/' $file_limma_in)

########## 5.2 submit jobs to run the scripts ########
#cd $shDir
# ls limma.$cellType.*.sh | awk '{print "sbatch "$0}' | bash
#cd $cDir
##################################################

######################## 6. limma to sce ########################
### please check that the format of data.id is: ^(cancerType).(dataset)$, and ther are no "." (dots) in cancerType and dataset 
### signature gene difference of BRCA.ElhamAzizi2018_InDrop
###    | sed 's/BRCA.ElhamAzizi2018.InDrop/BRCA.ElhamAzizi2018_InDrop/' \
join -1 1 -2 1 \
    <(cut -f 1-3 $file_int_in|sort -k 1r,1) \
    <(ls $cDir/OUT.int.$cellType.2nd/limma.sc/*/*.de.out.rda | perl -ane 'BEGIN{print "data.id\tdfile\n" } chomp; /.+limma.sc.(.+?).de.out/;print "$1\t$_\n"' | sort -k 1r,1) \
    | awk '!/(BRCA.ElhamAzizi2018.InDrop|HC.JiyuanZhang2020)/' \
    | sed 's/\s\+/\t/g' \
    | sed 's/BRCA.ElhamAzizi2018.10X/BRCA.ElhamAzizi2018_10X/' \
    > $file_limma_out

$sDir/wrapper.convertLimmaToSCE.R \
    --limmaFile $file_limma_out \
    --outPrefix $cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc \
    --ncores 8

######################## 7. prepare data for web ########################
mkdir $cDir/OUT.data.web
(
cat <<-HERE
$cDir/OUT.int.$cellType.2nd/int.$cellType.meta.tb.rds 
$cDir/OUT.int.$cellType.2nd/int.$cellType.sce.merged.rds
$cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.sce.pb.rds
$cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.gene.desc.tb.rds
$cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.geneTableLong.rds
$cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.geneTableLong.collapsed.rds 
$cDir/OUT.int.$cellType.2nd/int.$cellType.colSet.rds
HERE
) | perl -ane 'chomp; ($dname,$bname)=/^(.+)\/(.+)$/; print "ln -s $_ '$cDir'/OUT.data.web/'$g_prj_id'.$bname\n"' \
  | bash

#,DatasetName,DatasetSource,perMiniCluster,perMetaCluster,meta.perCell,geneTableLong,geneDesc,colSet
printf "$g_prj_id.$cellType,$g_prj_id.$cellType,$cellType,$g_prj_id.int.$cellType.sce.merged.rds,$g_prj_id.int.$cellType.limma.sc.sce.pb.rds,$g_prj_id.int.$cellType.meta.tb.rds,$g_prj_id.int.$cellType.limma.sc.geneTableLong.collapsed.rds,$g_prj_id.int.$cellType.limma.sc.gene.desc.tb.rds,$g_prj_id.int.$cellType.colSet.rds\n" > $cDir/OUT.data.web/dataset_map.csv

############### final results ###################
# $cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.gene.desc.tb.rds
# $cDir/OUT.int.$cellType.2nd/int.$cellType.limma.sc.sce.pb.rds
# $cDir/OUT.int.$cellType.2nd/int.$cellType.meta.tb.rds
# $cDir/OUT.int.$cellType.2nd/int.$cellType.seu.merged.rds
# $cDir/OUT.int.$cellType.2nd/int.$cellType.sce.merged.rds

### important columns: geneID median.F.rank
# OUT.int.$cellType.2nd/int.$cellType.gene.rank.tb.flt.rds



