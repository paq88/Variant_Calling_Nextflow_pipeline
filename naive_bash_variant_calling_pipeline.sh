#!/bin/bash
### BASIC PIPELINE ###
### NextFlow project ### 12.02.2023 ### Magda Mielczarek

### env/software
#conda activate nflow_env

### global parameters
t=5

### global vars
aln=/media/DANE/home/nflow/nextflow_project/data/ref/Bos_taurus.ARS-UCD1.2.dna.primary_assembly.25.fa
RG_file=/media/DANE/home/nflow/nextflow_project/data/reads_trimmed/NADIR_RG_2021.uniq.txt

in_dir=/media/DANE/home/nflow/nextflow_project/data/reads_trimmed/HOLPOLF005074063739
out_dir=/media/DANE/home/nflow/nextflow_project/outdir

echo "### Quality control"
time fastqc -t $t $in_dir/*.fastq.gz
#The output is located in input directory by default
#Can be changed by fastqc -t $t -o $out_dir_define_first $input

echo "### Alignment to the reference genome"
for fq1 in $in_dir/*trimm.forw.paired.fastq.gz
do
	
	bam_dir=$out_dir/$(dirname $fq1 | awk -F [/] '{print $9}')
	mkdir -p $bam_dir
	fq2=$(echo $fq1 | sed -e 's/forw/reve/')
	echo $fq1
	echo $fq2
	bam=$(basename $fq1 .trimm.forw.paired.fastq.gz).bwamem.bam
	tag=$(zcat $fq1 | head -n 1 | cut -d : -f 3,4 | sed 's/:/_/g')
	RG=$(grep $tag $RG_file)
	
	time bwa mem -t $t \
	-R $RG \
	$aln $fq1 $fq2 | \
	samtools view -bSh - > $bam_dir/$bam
done

echo "### Post-alignment"
# All bams belonging to 1 animal have to be merged.
mbam=$(basename $bam_dir)
time samtools merge -c -f $bam_dir/$mbam.merged $bam_dir/*.bwamem.bam
time samtools fixmate -m  $bam_dir/$mbam.merged $bam_dir/$(basename $mbam).fixmate.bam
time samtools sort -@0 -o $bam_dir/$(basename $mbam).fixmate.srt.bam $bam_dir/$(basename $mbam).fixmate.bam
time samtools markdup -r -s $bam_dir/$(basename $mbam).fixmate.srt.bam $bam_dir/$(basename $mbam).fixmate.srt.markdup.bam
time samtools index $bam_dir/$(basename $mbam).fixmate.srt.markdup.bam
time samtools flagstat $bam_dir/$(basename $mbam).fixmate.srt.markdup.bam > $bam_dir/$(basename $mbam).fixmate.srt.markdup.flagstat

echo "### Unnecessary files removal"
rm $bam_dir/*.bwamem.bam
rm $bam_dir/$mbam.merged
rm $bam_dir/$(basename $mbam).fixmate.bam
rm $bam_dir/$(basename $mbam).fixmate.srt.bam

echo "### SNP calling"
time bcftools mpileup -Ou -f $aln $bam_dir/$(basename $mbam).fixmate.srt.markdup.bam | bcftools call -vmO z -o $bam_dir/$(basename $mbam).fixmate.srt.markdup.vcf.gz

