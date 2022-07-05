#!/bin/bash

sDir=`dirname $0`
#iniFile="$sDir/../parameter/init_human.sh"
#_refData="/WPSnew/zhenglt/00.database/broad/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
optT=8

while getopts c:t: opt
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
		optT=$OPTARG
		;;
	'?')
		echo "Usage: $0 invalid option -$OPTARG" 
	    echo "Usage: $0 [-c iniFile] [-t threads, default 8] <sampleID> <outDir> <inBam> <inBarcode>"
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

if [ $# -lt 4 ]
then 
	echo "Usage: $0 [-c iniFile] [-t threads, default 8] <sampleID> <outDir> <inBam> <inBarcode>"
	exit 1
fi

echo begin at: `date`

#source $iniFile

sampleID=$1
outDir=$2
inBam=$3
inBarcode=$4

module load cellsnp/github

echo begin at: `date`

mkdir -p $outDir
mkdir -p $outDir/mquad
cd $outDir

cellsnp-lite -s $inBam -b $inBarcode --outDir=$outDir --chrom=chrM --UMItag Auto --minMAF 0 --minCOUNT 0 --genotype --gzip -p $optT

mquad --cellData $outDir --outDir $outDir/mquad -p $optT --minDP 5

$sDir/wrapper.vireoSNP.py -i $outDir/mquad -o $outDir/mquad/vireo

echo end at: `date`



