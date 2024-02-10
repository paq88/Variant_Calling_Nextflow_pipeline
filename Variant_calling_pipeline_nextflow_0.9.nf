#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//Parameters
// Here you should put your paths to out-directory, reference genome, RG file and reads
params.outdir = "/media/DANE/home/nflow/nextflow_project/programs/nextflow_test/out"
params.ref = "/media/DANE/home/nflow/nextflow_project/data/ref/Bos_taurus.ARS-UCD1.2.dna.primary_assembly.25.fa"
params.RG = "/media/DANE/home/nflow/nextflow_project/data/reads_trimmed/NADIR_RG_2021.uniq.txt"
params.reads = "/media/DANE/home/nflow/nextflow_project/data/reads_trimmed/HOLPOL*"

// You should if needed define number of cores in processes 

//Channels


ref_genome_ch = Channel.value( "${params.ref}")

RG_file_ch = Channel.value( "${params.RG}")

forw_reads_ch = Channel.fromPath("${params.reads}/*forw*fastq.gz", checkIfExists: true)

all_reads_ch = Channel.fromPath("${params.reads}/*fastq.gz", checkIfExists: true)


//Cosmetics
println """\
         V A R I A N T - C A L L I N G - N F - P I P E L I N E
         ===================================================
		 Authors: M.Sztuka, P.Hajduk, J.Liu
		 outdir       : ${params.outdir}
         """
		 .stripIndent()



//Proceses

process tuple_gen{
input: val(reads_path)
output:
tuple env(animal_id),env(full_path), emit: tuple_id_path

println"Tupple generation - forward reads"

script:
"""

animal_id=`echo ${reads_path} | cut -d "/" -f 9`
full_path=$reads_path


"""
}

process tuple_gen_qc{
input: val(reads_path)
output:
tuple env(animal_id),env(full_path), emit: tuple_id_path

println"Tupple generation - all reads"

script:
"""

animal_id=`echo ${reads_path} | cut -d "/" -f 9`
full_path=$reads_path


"""
}




process quality_controll{ 
cpus 5
maxForks 5
publishDir = "${params.outdir}/fastqc_out"
input:
tuple val(animal_id),path(reads)
output: path("*")

println"Quality control"

script:
"""
name=\$(basename $reads)
fastqc -c ${task.cpus} $reads > ${animal_id}_\${name}.html

"""
}



process read_alignment { 
publishDir = '/media/DANE/home/nflow/nextflow_project/programs/nextflow_test/out/bam_out'
cpus 5
maxForks 5

input:
tuple val(animal_id),val(forw)
val (ref)
val (RG_head)

output:

tuple val(animal_id), path('*.bam')    

println"Aligment"

script:
	"""



	reve=\$(echo $forw | sed -e 's/forw/reve/')


	#tag
	#wyciaganie ze sciezki nazwy pliku 
	tag=\$(basename -- "$forw")

	#usuwanie rozszerzen 
	tag="\${tag%%.*}"

	#modyfikacja - wyciaganie wylacznie 4 i 5 pola z tag
	tag=\$(echo \$tag | awk -F '_' '{print \$4,\$5}') 

	#zamiana spacji na _ 
	tag=\${tag/ /_}


	#plik RG
	head_tag=\$(zcat $forw | head -n 1 | cut -d : -f 3,4 | sed 's/:/_/g')

	RG=\$(grep \$head_tag $RG_head)
	#echo 'RG is'
	#echo \$RG

	bwa mem -t ${task.cpus} -R \$RG $ref $forw \$reve | samtools sort -n - > res_\${tag}.bam

"""
}


process post_alignment {

	maxForks 5
	cpus 5
	publishDir = '/media/DANE/home/nflow/nextflow_project/programs/nextflow_test/out/post_aln_out/indexed'


	input:
	val (bam_files)


	output:
   	path("*.markdup"), emit: markdup_samples

	println"Post_aligment"

	script:
	"""
	echo $bam_files

	
	mbam=\$(basename $bam_files)
	

	echo \$mbam
	samtools merge -@ ${task.cpus}-cn -f \$mbam.merged -b $bam_files 
	samtools fixmate -@ ${task.cpus} -m  \$mbam.merged \$mbam.fixmate 
	samtools sort -@ ${task.cpus} -o \$mbam.srt \$mbam.fixmate 
	samtools markdup -@ ${task.cpus} -r -s \$mbam.srt \$mbam.markdup 
	samtools index -@ ${task.cpus} \$mbam.markdup


	"""

}

process calling_SNP{
	publishDir = "${params.outdir}/VCF_out"
	cpus 5
	maxForks 5

	input:
	path(reads)
	path(ref_genome)
	output:
	path("*")

	script:
	"""
	name=\$(basename -- "$reads")
	bcftools mpileup -Ou -f $ref_genome $reads | bcftools call --threads ${task.cpus} -vmO z -o \$name.vcf.gz
	
	"""

}


//Workflow
workflow{

	//Tuples
id_tuple_forw=tuple_gen(forw_reads_ch).tuple_id_path
id_tuple_all=tuple_gen_qc(all_reads_ch).tuple_id_path

	//Quality control
quality_controll(id_tuple_all)

	//Aligment
bam_samples=read_alignment(id_tuple_forw,ref_genome_ch,RG_file_ch).collectFile(){item -> ["${item[0]}.txt","${item[1]}"+"\n"]}.view()


	//Post aligment
post_reads=post_alignment(bam_samples).markdup_samples

	//Calling SNP vcf generation
calling_SNP(post_reads,ref_genome_ch)

}

