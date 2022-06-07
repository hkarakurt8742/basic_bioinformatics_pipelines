#! /usr/bin/bash

sample_name=
fast5_loc=
guppy_out=
cpu_caller=
num_of_callers=
basecalling_config=
CPU=
ref_loc=

guppy_basecaller -i $fast5_loc -s $guppy_out --cpu_threads_per_caller $cpu_caller --num_callers $num_of_callers -c $basecalling_config

cat $guppy_out/\*.fastq > $guppy_out\/all_guppy.fastq

minimap2 -ax map-ont -t $CPU $ref_loc $guppy_out\/all_guppy.fastq > $sample_name\.sam

samtools view -@ $CPU $sample_name\.sam -o $sample_name\.bam

samtools sort -@ $CPU $sample_name\.bam -o $sample_name\_sorted.bam

samtools index -@ $CPU $sample_name\_sorted.bam

longshot --bam $sample_name\_sorted.bam --ref $ref_loc --out $sample_name\.vcf --snps


