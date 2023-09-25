#!/bin/bash

# split bam by CHR
outfold=$1;
seq_type=$2
input=$3


# output folder name
j2=$(basename $outfold)

# input file name
infile=$input"/"$j2"/"$j2".Split.bam"
echo ""
echo "input_file:" $infile
echo "output_directory:" $outfold

# Get the list of chromosomes for which we have read alignments
arr=($(samtools idxstats $infile | awk '{if($3 > 10) print $1}'))
echo "Splitting .bam files for following chromosomes: " ${arr[@]}
echo ""

# split the bam file chromosome wise
mkdir -p ${outfold}
 if [ $seq_type == "single" ]
   then
     for chrom in ${arr[@]}; do samtools view -b ${infile} ${chrom} > ${outfold}/$j2"."${chrom}.bam ; done;
   else
     for chrom in ${arr[@]}; do samtools view -b -f 0x2 ${infile} ${chrom} > ${outfold}/$j2"."${chrom}.bam ; done;
 fi

# create bam index
 for bam in ${outfold}/*.bam; do
         name=`basename $bam .bam`;
         samtools index $bam ${outfold}/${name}.bai; done;

