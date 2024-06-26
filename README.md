You can learn more in our article: 

[Nextflow vs. plain bash: different approaches to the parallelization of SNP calling from the whole genome sequence data](https://academic.oup.com/nargab/article/6/2/lqae040/7659592 )

Three pipelines were used for comparison 
- Native Bash application
- Single-process Nextflow where all of the pipeline was crammed into single Nextflow process
- Multi-process Nextflow where each step of the pipeline was seperate Nextflow process (the correct way to design pipeline in managed system such as Nextflow) "Variant_calling_pipeline_nextflow_0.9.nf"

Files used can be found on ncbi database here are coresponding accesion numbers with links:
fastq files for individuals:
- SRX2455298  https://www.ncbi.nlm.nih.gov/sra/SRX2455298[accn]
- SRX2455321  https://www.ncbi.nlm.nih.gov/sra/SRX2455321[accn]
- SRX2455312  https://www.ncbi.nlm.nih.gov/sra/SRX2455312[accn]
- SRX2455300  https://www.ncbi.nlm.nih.gov/sra/SRX2455300[accn]
- SRX2455319  https://www.ncbi.nlm.nih.gov/sra/SRX2455319[accn]


