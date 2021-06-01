#!/bin/bash

sDir=`dirname $0`
iniFile="$sDir/../parameter/init_human.sh"
_refData="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
###_refData="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optT="--localcores 4"
optM="--localmem 16"
optS="human"
optP=""

while getopts c:t:m:s:p: opt
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
	t)
		optT="--localcores $OPTARG"
		;;
	m)
		optM="--localmem $OPTARG"
		;;
    s)
        optS=$OPTARG
        ;;
    p)
        optP=$OPTARG
        ;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
		echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 16] [-s species, default human] [-p cellranger parameters] <outDir> <inDir> <sampleID>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 3 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 4] [-m memrory(GB), default 16] [-s species, default human] [-p cellranger parameters] <outDir> <inDir> <sampleID>"
	exit 1
fi

echo begin at: `date`

source $iniFile

outDir=$1
inDir=$2
sampleID=$3

if [ "$optS" == "human" ];then
    transcriptomeDir="/WPSnew/zhenglt/00.database/ensemble/10X/refdata-cellranger-GRCh38-1.2.0"
elif [ "$optS" == "mouse" ];then
    transcriptomeDir="/WPSnew/zhenglt/00.database/ensemble/10X/refdata-cellranger-mm10-1.2.0"
fi

mkdir -p $outDir
###export MODULESHOME=/usr/share/Modules
###. /usr/share/Modules/init/bash
###export MODULEPATH="/Share/BP/zhenglt/05.setting/modulefiles":/usr/share/Modules/modulefiles:/etc/modulefiles
#module load cellranger/2.1.1
module load cellranger/3.0.0

#outDir=`pwd`/OUT.cellranger
#mkdir -p $outDir
#mkdir -p $outDir/fq

echo begin at: `date`
cd $outDir

echo cellranger count --id=$sampleID $optT $optM \
    --fastqs=$inDir \
    --sample=$sampleID \
    --transcriptome=$transcriptomeDir $optP


cellranger count --id=$sampleID $optT $optM \
    --fastqs=$inDir \
    --sample=$sampleID \
    --transcriptome=$transcriptomeDir $optP

    ####--expect-cells=3000
echo end at: `date`



