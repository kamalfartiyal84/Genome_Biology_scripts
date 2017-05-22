#!/bin/bash

bin_size=$1
sample_sheet=$2
bam_dir=$3
aligner=$4
mem=$5

if [[ -z $bin_size || -z $sample_sheet || -z $bam_dir ]]
then
	echo "Usage: `basename $0` bin_size sample_sheet bam_dir [aligner] [mem]" 1>&2
	exit 1
fi 

if [[ -z $mem ]]
then
	mem=4096
fi

if [[ -z $aligner ]]
then
	aligner=""
else
	aligner=".${aligner}"
fi

bam_dir=`echo $bam_dir | sed "s|/$||"`

if [ ! -d $bam_dir ]
then
	echo "Directory $bam_dir not found"
	exit 1
fi

mkdir -p logs

IFS="	"
sed 1d $sample_sheet | while read id sample
do
	if [ ! -e ${id}.binSize${bin_size}.copyNumber.txt ]
	then
		echo $id

		bsub -q bioinformatics -R "rusage[mem=${mem}]" -o logs/${id}.binSize${bin_size}.qdnaseq.%J.log -J qdnaseq_${id} << EOF
/home/software/R/R-3.2.3/bin/Rscript `dirname $0`/qdnaseq.R ${bam_dir}/${id}${aligner}.bam ${id} "${sample}" ${id}.binSize${bin_size} $bin_size
EOF

	fi
done

