'''
RNA_seq_STAR_RSEM.snakefile
Kiran Narta

 Map paired end RNA seq reads with STAR and quantify transcripts with RSEM.
--------------------------------------------------------------------------
#### inspired from : https://github.com/slowkow/snakefiles/blob/master/star_express/Snakefile

Requirements: STAR, RSEM, sambamba, multiqc

Usage: snakemake --latency-wait 60 --cluster-config cluster.json --jobs 100 --cluster "sbatch -t {cluster.time} -p {cluster.partition} -n {cluster.n} --mem={cluster.mem} --output=out_err/%j.out --error=out_err/%j.err --cpus-per-task={cluster.cpus-per-task}" -s Snakefile
 
'''

import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
from glob import glob

#### Globals --------------------------------------------------------
STARGenome = "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/STAR_Genome" ## path to indexed STAR genome
RSEMGenome = "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/RSEM_Genome"	## path to RSEM genome
GTF = "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/gencode.vM10.annotation_plV430GFeatures1.gtf" ## path to config file
out_dir = "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/SyngeneicModels/Batch1/"    ## path to all out files

FILES = json.load(open("samples.json"))
SAMPLES = sorted(FILES.keys())

###### Function to remove suffix from a string ###
def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

##### RULES ---------------------------------------------------------------

rule all:
	input:
		join(out_dir, 'SJ.out.pass1_merged.tab'),
		expand(join(out_dir, '{sample}', 'RSEM', '{sample}_rsem.genes.results'), sample=SAMPLES)


# Merge multiple FASTQ files that correspond to a single sample.
# Also, ensure that the merged FASTQ has correctly paired mates.
# Unpaired mates are excluded from the merged FASTQ.
rule merge_fastqs:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        r1 = temp(join(out_dir, '{sample}', '{sample}_R1.fastq.gz')),
        r2 = temp(join(out_dir, '{sample}', '{sample}_R2.fastq.gz'))
    run:
        import gzip
        # Open the output file handles.
        out1 = gzip.open(output.r1, 'wb')
        out2 = gzip.open(output.r2, 'wb')
        # Yield a tuple for each read in the FASTQ file.
        def fastq_records(filename):
            with gzip.open(filename) as f:
                while True:
                    try:
                        name = f.readline()
                    except:
                        break
                    seq = f.readline()
                    sep = f.readline()
                    qual = f.readline()
                    yield (name.split()[0], name + seq + sep + qual)
        # Loop through all (R1,R2) pairs.
        for i in range(len(input.r1)):
            in1 = fastq_records(input.r1[i])
            in2 = fastq_records(input.r2[i])
            d1 = {}
            d2 = {}
            done1 = False
            done2 = False
            while not done1 or not done2:
                try:
                    n1, r1 = next(in1)
                    d1[n1] = r1
                except:
                    done1 = True
                try:
                    n2, r2 = next(in2)
                    d2[n2] = r2
                except:
                    done2 = True
                if not done1 and n1 in d2:
                    out1.write(d1[n1])
                    out2.write(d2[n1])
                    d1.pop(n1)
                    d2.pop(n1)
                if not done2 and n2 in d1:
                    out1.write(d1[n2])
                    out2.write(d2[n2])
                    d1.pop(n2)
                    d2.pop(n2)
        # Close the output file handles.
        out1.close()
        out2.close()



## STAR-Pass1: Map paired end reads to genome, Delete output bam file
rule star_pass1:
	input:
		r1 = rules.merge_fastqs.output.r1,
		r2 = rules.merge_fastqs.output.r2,
		genomedir = STARGenome,
		gtf = GTF
	output:
		bam=temp(join(out_dir, '{sample}', 'pass1','Aligned.out.bam')),
		sj= join(out_dir, '{sample}','pass1', 'SJ.out.tab')
	log:
		join(out_dir, '{sample}', 'pass1', 'star.map.log')
	benchmark:
		join(out_dir, '{sample}', 'pass1', 'star.map.benchmark.tsv')
	threads:
		10
	resources:
		mem=40000
	run:
		shell('cd ' + join(out_dir, '{wildcards.sample}', 'pass1') +
			'&& STAR'
			' --runThreadN {threads}'
			' --genomeDir {input.genomedir}'
			' --sjdbGTFfile {input.gtf}'
			' --outSAMtype BAM Unsorted'
			' --outSAMunmapped Within'
			' --readFilesCommand zcat'
			' --readFilesIn {input.r1} {input.r2}'
			' --outFileNamePrefix ' + join(out_dir, '{wildcards.sample}', 'pass1') + '/'
			' > {log} 2>&1')

## Merge Splice junctions from all files
rule merge_splice_junctions:
	input:
        	sjs = expand(join(out_dir, '{sample}', 'pass1', 'SJ.out.tab'), sample = SAMPLES)
	output:
		sjs = join(out_dir, 'SJ.out.pass1_merged.tab')
	log:
		join(out_dir, 'merge_splice_junctions.log')
	shell:
		# Retain splice junctions with at least 3 uniquely mapped fragments	# per sample.
		"cat {input.sjs} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output.sjs}"


## Generate Genome
rule GenomeGenerate:
	input:
		genomefasta = STARGenome,	
		gtf = GTF,
		sjs = rules.merge_splice_junctions.output.sjs
	output:
		index = join(dirname(GTF), 'STARGenomePass2', 'Genome')
	log:
		join(dirname(GTF), 'STARGenomePass2', 'star.index.log')
	benchmark:
		join(dirname(GTF), 'STARGenomePass2', 'star.index.benchmark.tsv')
	threads:
		10
	resources:
		mem= 40000
	run:
		shell('STAR'
			' --runMode genomeGenerate'
			' --genomeDir ' + join(dirname(GTF), 'STARGenomePass2') +
			' --genomeFastaFiles {input.genomefasta}'
			' --sjdbGTFfile {input.gtf}'
			' --sjdbOverhang 100'
			' --runThreadN {threads}'
			' --sjdbFileChrStartEnd {input.sjs}'
			' > {log} 2>&1')
			


### STAR-Pass2: Map PE RNAseq reads, make bam relative to transcripts and genome, Count reads mapping to gene, count reads supporting each splice junction
rule star_pass2:
	input:
		r1 = rules.merge_fastqs.output.r1,
                r2 = rules.merge_fastqs.output.r2,
                genomedir = rules.GenomeGenerate.output.index
	output:
		bam=join(out_dir, '{sample}', 'pass2','Aligned.out.bam'),
		t_bam=join(out_dir,'{sample}', 'pass2', 'Aligned.toTranscriptome.out.bam'),
		counts=join(out_dir, '{sample}', 'pass2', 'ReadsPerGene.out.tab'),
		sj= join(out_dir, '{sample}', 'pass2', 'SJ.out.tab')
	log:
		join(out_dir, '{sample}', 'pass2', 'star.map.log')
	benchmark:
		join(out_dir, '{sample}', 'pass2', 'star.map.benchmark.tsv')
	threads:
		10
	resources:
		mem=40000
	run:
		shell('cd ' + join(out_dir, '{wildcards.sample}', 'pass2')+
		'&& STAR'
		' --runThreadN {threads}'
		' --genomeDir ' + dirname(input.genomedir) +
		' --readFilesCommand zcat'
		' --readFilesIn {input.r1} {input.r2}'
		' --quantMode TranscriptomeSAM GeneCounts'
		' --outSAMtype BAM Unsorted'
		' --outSAMunmapped Within'
		' --outFileNamePrefix ' + join(out_dir, '{wildcards.sample}', 'pass2') + '/'
		' --outFilterMultimapNmax 1'
		' > {log} 2>&1')


## Sambamba Sort
rule sort_bam:
	input:
                tbam = rules.star_pass2.output.t_bam
        output:
                ts_bam = join(out_dir,'{sample}', 'pass2', 'Aligned.toTranscriptome.out.sorted.bam')
        threads:
                10
        run:
                shell('cd ' + join(out_dir, '{wildcards.sample}', 'pass2')+
                "&& sambamba sort -p -t {threads} -o {output.ts_bam} {input.tbam}"
		'> {log} 2>&1')

#Convert bam to RSEM bam
rule converttorsem:
	input:
		sorted_tbam = rules.sort_bam.output.ts_bam
	output:
		rsem= join(out_dir, '{sample}', 'RSEM', '{sample}_Aligned.toTranscriptome.out.sorted_RSEM.bam')
	threads:
		10
	run:
		shell('cd ' + join(out_dir, '{wildcards.sample}', 'RSEM')+
		"&& convert-sam-for-rsem -p {threads} {input.sorted_tbam} {output.rsem}"
		'> {log} 2>&1')


# Run RSEM  
rule run_rsem:
	input:
		rsembam = rules.star_pass2.output.t_bam,
		genome= multiext("/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/RSEM_Genome/GRCm10_plV430G", ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa")
		#genome= "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/RSEM_Genome/GRCm10_plV430G/GRCm10_plV430G" 
	output:
		genes = join(out_dir, '{sample}', 'RSEM', '{sample}_rsem.genes.results'),
		isoforms = join(out_dir, '{sample}', 'RSEM', '{sample}_rsem.isoform.results'),
		statcnt= join(out_dir, '{sample}', 'RSEM', '{sample}.stat', '{sample}_rsem.cnt'),
		statmodel= join(out_dir, '{sample}', 'RSEM', '{sample}.stat', '{sample}_rsem.model'),
		stattheta= join(out_dir, '{sample}', 'RSEM', '{sample}.stat', '{sample}_rsem.theta')
	threads:
		10
	log:
		join(out_dir, '{sample}', 'RSEM', 'rsem.log')
	shell:
		"rsem-calculate-expression --bam --no-bam-output --append-names -p 10 --paired-end --forward-prob 0 {input.rsembam} /work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/RSEM_Genome/GRCm10_plV430G {out_dir}/{wildcards.sample}/RSEM/{wildcards.sample}_rsem >& {log}"

## te-expression --bam --no-bam-output --append-names -p 10 --paired-end --forward-prob 0 G287M0/pass2/Aligned.toTranscriptome.out.bam ../../Genome_mm10_plV430G/RSEM_Genome/GRCm10_plV430G G287M0/rsemtest/G287M0_rsem >& G287M0/rsemtest.log"Sort Bams



## Convert bam to Cram
#rule bamTOcram:
#	input:
#		starGenbam = rules.converttorsem.output.rsem,
#		genome = "/work/morrissy_lab/Kiran/Collabrative_Work/MikeMonument_MouseSarcoma/Genome_mm10_plV430G/STAR_Genome/Genome" 
#	output:
#	
		
		


