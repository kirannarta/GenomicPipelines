##### Snakemake pipeline for Variation calling Whole genome

# https://h3abionet.github.io/H3ABionet-SOPs/Variant-Calling#the-gatk-resource-bundle
#### GATK best practices: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

### Usage: snakemake --snakefile Germline_variant_calling_pipeline.snakefile --latency-wait 60 --cluster-config cluster.json --jobs 100 --cluster "sbatch -t {cluster.time} -p {cluster.partition} -n {cluster.n} --mem={cluster.mem} --output=out_err/%j.out --error=out_err/%j.out" 
#snakemake --snakefile Germline_variant_calling_pipeline.snakefile --latency-wait 5000 --cluster-config cluster.json --jobs 400 --cluster "sbatch -t {cluster.time} -p {cluster.partition} -n {cluster.n} --mem={cluster.mem}"

from os.path import join, dirname
from subprocess import check_output

workdir: "/work/morrissy_lab/Kiran/LIVIA_project/Medullo_PDXcellLine_WGS/processedData/Results_Med411/"

(SAMPLES,)=glob_wildcards("/work/morrissy_lab/collaborations/medullos_toronto/Mailbox/Med411/{sample}_R1.fastq.gz")

#Rules-----------------------
rule all:
    input:
        expand("Aligned/{sample}.hg38.bam", sample=SAMPLES),
	expand("Aligned/{sample}.sorted.hg38.bam", sample=SAMPLES),
	expand("Aligned/{sample}.sorted.hg38.bam.bai", sample=SAMPLES),
	expand("GATK/{sample}.sorted.dedup.hg38.indelrealined_BQSR.bam", sample=SAMPLES)


### FastQC

### BWA-MEM alignment https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html
rule alignment:
	input:
		bwa_genome="/bulk/morrissy_bulk/REFERENCES/human/hg38/genomes/BWAIndex/genome.fa",
		r1="/work/morrissy_lab/collaborations/medullos_toronto/Mailbox/Med411/{sample}_R1.fastq.gz",
		r2="/work/morrissy_lab/collaborations/medullos_toronto/Mailbox/Med411/{sample}_R2.fastq.gz"
	output:
		bam="Aligned/{sample}.hg38.bam"
	threads:10
	log:
		"logs/bwa_mem/{sample}.log"
	params:rg=r"@RG\tID:{sample}\tSM:{sample}"
	shell:
		"bwa mem -t {threads} -R '{params.rg}' {input.bwa_genome} {input.r1} {input.r2} > {output.bam}"		

### Picard Sort files
rule sort_bam:
	input:
		"Aligned/{sample}.hg38.bam"
	output:
		sorted_bam="Aligned/{sample}.sorted.hg38.bam"
	threads: 10
	shell:
		"samtools sort -@ {threads} {input} -o {output.sorted_bam}"
	
### Index bam files
rule index_bam:
	input:
		"Aligned/{sample}.sorted.hg38.bam"
	output:
		"Aligned/{sample}.sorted.hg38.bam.bai"
	shell:
		"samtools index {input}"

### Add read group ## GATK BQSR requires it
rule RG_sort:
        input:
                bam = rules.sort_bam.output.sorted_bam,
        output:
                RG = "Aligned/{sample}.sorted.RG.hg38.bam"
        shell:
                "picard AddOrReplaceReadGroups I= {input.bam} O= {output.RG} RGID={wildcards.sample} RGLB=WGS_DNA RGPL=illumina RGPU=machine RGSM={wildcards.sample}"


### Picard Mark Duplicates
rule Mark_duplicates:
	input:
		rules.RG_sort.output.RG
	output:
		dupmarked="Aligned/{sample}.sorted.RG.dedup.hg38.bam",
		metrics = "Aligned/{sample}.sorted.hg38.dedup.metrics"
	shell:
		"picard -Xmx8g MarkDuplicates I= {input} O= {output.dupmarked} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.metrics}"

### GATK: Indel Realignment
rule create_target:
        input:
                bam = rules.Mark_duplicates.output.dupmarked,
                genome = "/work/morrissy_lab/Kiran/RESOURCES/Refrences/Grch38/genome.fa",
                knownindel = "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        output:
                int = "GATK/{sample}.realign.intervals"
        shell:
                "GenomeAnalysisTK -Xmx4G -T RealignerTargetCreator -R {input.genome} -known {input.knownindel} -I {input.bam} -o {output.int}"


rule indel_realign:
        input:
                bam = rules.Mark_duplicates.output.dupmarked,
                genome = "/work/morrissy_lab/Kiran/RESOURCES/Refrences/Grch38/genome.fa",
                knownindel = "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
		int = rules.create_target.output.int
        output:
                bam = "GATK/{sample}.sorted.dedup.hg38.indelrealined.bam"
        shell:
                "gatk -Xmx4G -T IndelRealigner -R {input.genome} -known {input.knownindel} -I {input.bam} -targetIntervals {input.int} -o {output.bam}"

### GATK: BAse Recalibration
rule BQSR1:
        input:
                bam = "GATK/{sample}.sorted.dedup.hg38.indelrealined.bam",
                genome = "/work/morrissy_lab/Kiran/RESOURCES/Refrences/Grch38/genome.fa",
                knownsnp = "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                knownindels = "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz",
		knownsnp1 = "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf",
		knownindels2= "/work/morrissy_lab/Kiran/RESOURCES/GATK_resources/hg38_resources/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

        output:
                table = "GATK/{sample}.recal_data.table"
        shell:
                "gatk -Xmx4G -T BaseRecalibrator -R {input.genome} -I {input.bam} -knownSites {input.knownsnp} -knownSites {input.knownindels} -knownSites {input.knownsnp1} -knownSites {input.knownindels2} -o {output.table}"

rule BQSR2:
        input:
                genome = "/work/morrissy_lab/Kiran/RESOURCES/Refrences/Grch38/genome.fa",
                table = rules.BQSR1.output.table,
                bam = rules.indel_realign.output.bam
        output:
                bam = "GATK/{sample}.sorted.dedup.hg38.indelrealined_BQSR.bam"
        shell:
                "gatk -Xmx4G -T PrintReads -R {input.genome} -I {input.bam} -BQSR {input.table} -o {output.bam}"


#### Mutect2 pairwise 

### Variant Annotation
#run in command line:=>  while read sample; do sbatch -n 2 --mem 40000 --time=5-24:00:00 --partition=cpu2019 -o out_err/annovar.%j.out -e out_err/annovar.%j.err --wrap="perl /work/morrissy_lab/Kiran/RESOURCES/annovar/table_annovar.pl ../${sample}.sorted.chr17.dedup.hg38.indelrealined_postrecal_BQSR_raw.vcf /work/morrissy_lab/Kiran/RESOURCES/annovar/humandb/ -buildver hg38 -out ${sample}.chr17_annotation -remove -otherinfo -protocol refGene,gnomad211_exome,avsnp150,cosmic70,icgc28,ljb26_all,clinvar_20210123 -operation g,f,f,f,f,f,f -nastring . -vcfinput" ; done < vcf_files.txt  

