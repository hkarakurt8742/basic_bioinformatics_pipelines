#! /usr/bin/bash

file_loc=
sra_name=
genome_dir=
read_info=
CPU=
known_sites=
bed_location=
gatk_loc=
annovar_db_loc=
annovar_loc=
sra_toolkit_loc=

$sra_toolkit_loc\/fasterq-dump --threads $CPU --progress --split-3 $sra_name

bwa mem -M -t $CPU -R $read_info $genome_dir $file_loc/$sra_name\_1.fastq $file_loc/$sra_name\_2.fastq  | samtools sort -@ $CPU -o $sra_name\.bam -

samtools index -@ $CPU $sra_name\.bam

rm $file_loc/$sra_name\_1.fastq
rm $file_loc/$sra_name\_2.fastq

picard MarkDuplicates INPUT=$sra_name\.bam OUTPUT=$sra_name\.marked.bam METRICS_FILE=$sra_name\_metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

$gatk_loc BaseRecalibrator -I $sra_name\.marked.bam  -R $genome_dir --known-sites $known_sites -O $sra_name\_recal_data.table

$gatk_loc ApplyBQSR -R $genome_dir -I $sra_name\.marked.bam --bqsr-recal-file $sra_name\_recal_data.table -O $sra_name\.marked.recal.bam

$gatk_loc --java-options "-Xmx4g" HaplotypeCaller -R $genome_dir -I $sra_name\.marked.recal.bam -O $sra_name\.vcf.gz --dbsnp $known_sites -L $bed_location --native-pair-hmm-threads $CPU

$gatk_loc SelectVariants -R $genome_dir -V $sra_name\.vcf.gz --select-type-to-include SNP -O $sra_name\_snps.vcf.gz

$gatk_loc VariantFiltration -R $genome_dir -V $sra_name\_snps.vcf.gz -O $filtered_snps.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNPFilter"

$annovar_loc\table_annovar.pl $sra_name\_filtered_snps.vcf.gz $annovar_db_loc -buildver hg19 -out $sra_name\_filtered_snps_anno -remove -protocol genomicSuperDups,cpgIslandExt -operation r,r -arg "-colsWanted 27","-colsWanted 10" -nastring . -vcfinput -polish

