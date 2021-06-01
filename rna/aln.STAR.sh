#!/bin/bash -eu


sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"

optM=30
optMode=""

while getopts c:m: opt
do
	case $opt in 
	c)	
		if [ -f $OPTARG ]
		then
			iniFile="$OPTARG"
		else
			echo "WARNING: invalid reference file ($OPTARG), default will be used"
		fi
		;;
	m)
		optMode=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-m mode, default ''] <sampleID> <fq1> <fq2> <outDir>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-m mode, default ''] <sampleID> <fq1> <fq2> <outDir>"
	exit 1
fi

source $iniFile

STAR_GenomeDir=/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/star
gene_model_file=/WPSnew/zhenglt/00.database/gencode/v24/gencode.v24lift37.annotation.gtf
STAR_FUSION_DB_Dir=/WPSnew/zhenglt/00.database/tools/star-fusion/GRCh37_v24_mybuild/ctat_genome_lib_build_dir
module load STAR/2.6.1d
###module unload java/1.7.0_79
module load java/1.8.0_171
###module unload gatk/3.3-0
####module load gatk/3.8-0
module load subread/1.6.3
module load trimmomatic/0.33
module load blast/2.8.1+
module load STAR-Fusion/1.5.0

### with ERCC
#####export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/withERCC.v1/hg19.order.after.gsnap.fa
### without ERCC
#####export REF=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/gmap-gsnap/hg19.order.after.gsnap.fa
#####knownSites=/DBS/DB_temp/zhangLab/broad/bundle/2.8/hg19/dbsnp_138.hg19.reorder.vcf

sampleID=$1
fq1=$2
fq2=$3
outDir=$4
####REF=$5

### for most analysis, parallel mode is not supported currently 
optNT=1
###optNT=8

optL=''
optM=`echo "scale=0;$optM/1.5" | bc`
max_reads=`echo 250000*$optM | bc`

mkdir -p $outDir

newFQ1=$fq1
newFQ2=$fq2
len_trim_to=101

if [ "$optMode" == "trim" ];then
    if [ ! -f "$outDir/$sampleID.clean.P.R1.fq.gz" ]
    then
    	java -Xmx${optM}g -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar PE \
    		-threads 8 \
    		$fq1 \
    		$fq2 \
    		$outDir/$sampleID.clean.P.R1.fq.gz \
    		$outDir/$sampleID.clean.UP.R1.fq.gz \
    		$outDir/$sampleID.clean.P.R2.fq.gz \
    		$outDir/$sampleID.clean.UP.R2.fq.gz \
    		ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/NexteraPE-PE.fa:2:30:10 \
    		CROP:$len_trim_to \
    		TRAILING:3 \
    		MAXINFO:50:0.25 \
    		MINLEN:36 \
    		TOPHRED33 
    	echo "trimmomatic done. (" `date` ")"
    fi
    newFQ1=$outDir/$sampleID.clean.P.R1.fq.gz
    newFQ2=$outDir/$sampleID.clean.P.R2.fq.gz
fi

STAR \
    --runThreadN 8 \
    --genomeDir $STAR_GenomeDir \
    --readFilesIn $newFQ1 $newFQ2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $outDir/$sampleID. \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outSAMattrRGline ID:$sampleID CN:BIOPIC LB:$sampleID PL:illumina PU:$sampleID SM:$sampleID \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped None \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 12 \
    --alignSJDBoverhangMin 10 \
    --alignMatesGapMax 100000 \
    --alignIntronMax 100000 \
    --chimSegmentReadGapMax 3 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --outSAMstrandField intronMotif \
    --chimOutJunctionFormat 1

samtools index $outDir/${sampleID}.Aligned.sortedByCoord.out.bam

echo "... using $max_reads reads in memory (parameter optM: $optM*1.5)"
java -Xmx${optM}g -jar $PICARD/picard.jar SortSam \
		I=$outDir/${sampleID}.Aligned.sortedByCoord.out.bam \
		O=$outDir/${sampleID}.Aligned.sortedByRName.out.bam \
		MAX_RECORDS_IN_RAM=$max_reads \
		TMP_DIR=$outDir \
		SO=queryname \
		VALIDATION_STRINGENCY=SILENT
###samtools index $outDir/${sampleID}.Aligned.sortedByRName.out.bam

####### gene expression
featureCounts -p -T 2 \
    -a $gene_model_file \
    -o $outDir/$sampleID.subread.exp \
    $outDir/${sampleID}.Aligned.sortedByRName.out.bam
    ##$outDir/${sampleID}.Aligned.sortedByCoord.out.bam

featureCounts -p -T 2 -O \
    -a $gene_model_file \
    -o $outDir/$sampleID.subread.exp.optO \
    $outDir/${sampleID}.Aligned.sortedByRName.out.bam
    ##$outDir/${sampleID}.Aligned.sortedByCoord.out.bam

rm $outDir/${sampleID}.Aligned.sortedByRName.out.bam    

####### fusion
STAR-Fusion --genome_lib_dir $STAR_FUSION_DB_Dir \
    -J $outDir/$sampleID.Chimeric.out.junction \
    --CPU 8 \
    --output_dir $outDir/star-fusion

exit

####### variant calling ######

inBam=$outDir/${sampleID}.Aligned.sortedByCoord.out.bam
#####inBam=$outDir/${sampleID}.Aligned.sortedByCoord.out.RG.bam
outBam=$outDir/$sampleID.GATK.RNA.ready.bam

if [ -f $outBam ]
then
    echo "## $outBam exists, skip this step"
	exit 0
fi

if [ ! -f $inBam.bai ]
then
	samtools index $inBam
fi

echo ">>> Marking duplicates"
this_bam=$outDir/$sampleID.rmDup.bam
if [ ! -f $this_bam ]
then
    echo "using $max_reads reads in memory (parameter optM: $optM*1.5)"
    java -Xmx${optM}g -jar $PICARD/picard.jar MarkDuplicates \
    	TMP_DIR=$outDir \
    	I=$inBam \
    	O=$this_bam \
    	M=${this_bam%.bam}.metrics \
    	VALIDATION_STRINGENCY=SILENT \
    	ASSUME_SORTED=true \
    	REMOVE_DUPLICATES=false \
    	MAX_RECORDS_IN_RAM=$max_reads
    
    rm -f ${this_bam%.bam}.bai
    samtools index $this_bam
fi

echo ">>> splitNTrim"
optM=`echo "scale=0;$optM/1.5" | bc`
inBam=$this_bam
this_bam=$outDir/$sampleID.GATK.RNA.splitNCigar.bam
if [ ! -f $this_bam.bai ]
then
    java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
        -T SplitNCigarReads \
        -R $REF \
        -I $inBam \
        -o $this_bam \
        -U ALLOW_N_CIGAR_READS \
        -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
        ##-fixNDN
    ## for gsnap result, not nessesary:
    ##-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
    ## for SplitNCigarReads, parallel mode is not supported currently 
	##    -nt $optNT \
    samtools index $this_bam
fi

rm $inBam
rm $inBam.bai


this_bam=$outDir/$sampleID.GATK.RNA.realn.bam
if [ ! -f $this_bam.bai ]
then
    echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
    java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    	-T RealignerTargetCreator \
	    -I $outDir/$sampleID.GATK.RNA.splitNCigar.bam \
	    -R $REF \
	    -o $outDir/$sampleID.GATK.RNA.realn.intervals $optL \
	    -nt $optNT \
	    -rf BadCigar
    echo ">>> Running the realigner over the targeted intervals"
    java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
	    -T IndelRealigner \
	    -I $outDir/$sampleID.GATK.RNA.splitNCigar.bam \
	    -R $REF \
	    -o $outDir/$sampleID.GATK.RNA.realn.bam \
	    -targetIntervals $outDir/$sampleID.GATK.RNA.realn.intervals \
	    -LOD 5 $optL \
	    -nt $optNT \
	    -rf BadCigar
    samtools index $outDir/$sampleID.GATK.RNA.realn.bam
fi

echo $knownSites


this_bam=$outDir/$sampleID.GATK.RNA.recal.bam
if [ ! -f $this_bam.bai ]
then
    echo ">>> Counting covariates"
    optM=`echo "scale=0;$optM/1.2" | bc`
    java -Djava.io.tmpdir=$outDir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
    	-T BaseRecalibrator \
    	-R $REF \
    	--knownSites $knownSites \
    	--disable_indel_quals \
    	-I $outDir/$sampleID.GATK.RNA.realn.bam \
    	-o $outDir/$sampleID.GATK.RNA.recal.grp \
    	-cov ReadGroupCovariate \
    	-cov QualityScoreCovariate \
    	-cov CycleCovariate \
    	-cov ContextCovariate \
    	-rf BadCigar \
        --validation_strictness SILENT
    
    echo ">> Table recalibration"
    if [ "`grep -v '#' $outDir/$sampleID.GATK.RNA.recal.grp | grep -v "EOF" | wc -l`" = "1" ]
    then
        echo "no recal.grp"
        ln $outDir/$sampleID.GATK.RNA.realn.bam $outDir/$sampleID.GATK.RNA.recal.bam
    else
    	java -Djava.io.tmpdir=$outDir -Xms${optM}g -Xmx${optM}g -jar $GATK/GenomeAnalysisTK.jar \
    	-T PrintReads \
    	-R $REF \
    	-I $outDir/$sampleID.GATK.RNA.realn.bam \
    	-o $outDir/$sampleID.GATK.RNA.recal.bam \
    	-BQSR $outDir/$sampleID.GATK.RNA.recal.grp \
    	--disable_indel_quals \
    	-EOQ \
    	-rf BadCigar
    fi
    samtools index $outDir/$sampleID.GATK.RNA.recal.bam
fi

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $REF \
    -I $outDir/$sampleID.GATK.RNA.recal.bam \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
	-nt $optNT \
    --num_cpu_threads_per_data_thread 4 \
    -o $outDir/$sampleID.GATK.RNA.recal.var.vcf

java -Xms${optM}g -Xmx${optM}g -Djava.io.tmpdir=$outDir -jar $GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REF \
    -V $outDir/$sampleID.GATK.RNA.recal.var.vcf \
    -window 35 -cluster 3 \
    -filterName FS -filter "FS > 30.0" \
    -filterName QD -filter "QD < 2.0" \
    -o $outDir/$sampleID.GATK.RNA.recal.var.flt.vcf


echo end at: `date` 
