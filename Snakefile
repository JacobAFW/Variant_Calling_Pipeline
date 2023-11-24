# Snakefile 

# Prerequisites (on Gadi): source /g/data/pq84/malaria/snakemake_pipeline/snakemake_setup.sh

# Define paths for ref genomes using config file
configfile: "config/config.yaml"
fasta_path=config['fasta_path']
picard_path=config['picard_path']
gatk_path=config['gatk_path']
known_indels_path=config['known_indels_path']
known_sites_path=config['known_sites_path']
bed_file_path=config['bed_file_path']
source_dir=config['source_dir']

# Define input files - 

## Defining the samples to be used for the {sample} wildcard

SAMPLES, = glob_wildcards(f'{source_dir}/{{sample}}_1.fastq.gz') 

## Chromosome names for bed files and subsetting jobs
with open(config['bed_file_path']) as f:
    CHROMOSOME = f.read().splitlines()
    CHROMOSOME = [p.split('\t')[0] for p in CHROMOSOME]

# Define final files

rule all:
    input:
        expand("output/mapped_reads/{sample}.bam", sample = SAMPLES),
        expand("output/sorted/{sample}_sorted.bam", sample = SAMPLES),
        expand("output/sorted/{sample}_sorted.bam.bai", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked.bam", sample = SAMPLES),
        expand("output/bam_recal/{sample}_picard_metrics_file.txt", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked_reheader.bam", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked_reheader.bam.bai", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked_realigner.intervals", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked_realigned.bam", sample = SAMPLES),
        expand("output/bam_recal/{sample}_dupmarked_realigned_recal.table", sample = SAMPLES),
        expand("output/bam_recal/{sample}_recalibrated.bam", sample = SAMPLES),
        expand("output/calling/{sample}.g.vcf.gz", sample = SAMPLES),
        "output/calling/GATK_combined.g.vcf.gz",
        expand("data/chromosomes/{chromosome}.bed", chromosome = CHROMOSOME)

# Define local rules - not run with scheduler
localrules: all, bed_for_chrom

# Create bed for each chromosome/contig - to be used for joint calling

rule bed_for_chrom:
    input:
        bed=bed_file_path
    output:
        "data/chromosomes/{chromosome}.bed"
    run:
        with open(input.bed) as f:
            chrom_bed = f.read().splitlines()
            chrom_bed = [p.split('\t') for p in chrom_bed]

        for i in chrom_bed:
            content = '\t'.join(i)
            f = open(output[0], 'w')
            f.write(content)
            f.close()


# Combine gVCFs

rule combine_gvcfs:
    input:
        fasta=fasta_path
    output:
        "output/calling/GATK_combined.g.vcf.gz"
    params:
        gatk=gatk_path,
        gvcfs = lambda w: " -V " + " -V ".join(expand("output/calling/{sample}.g.vcf.gz", sample = SAMPLES))
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T CombineGVCFs \
        -R {input.fasta} \
        {params.gvcfs} \
        -o {output}
        """

# Call haplotypes - GATK

rule haplotype_caller:
    input:
        bam="output/bam_recal/{sample}_recalibrated.bam",
        fasta=fasta_path
    output:
        "output/calling/{sample}.g.vcf.gz"
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T HaplotypeCaller \
        -ERC GVCF \
        --minPruning 3 \
        --maxNumHaplotypesInPopulation 200 \
        --max_alternate_alleles 3 \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -contamination 0.0 \
        -G Standard \
        -R {input.fasta} \
        -I {input.bam} \
        -o {output}
        """

# Get recalibrated bams

rule print_reads:
    input:
        bam="output/bam_recal/{sample}_dupmarked_realigned.bam",
        table="output/bam_recal/{sample}_dupmarked_realigned_recal.table",
        fasta=fasta_path,
        bed=bed_file_path
    output:
        "output/bam_recal/{sample}_recalibrated.bam"
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T PrintReads \
        -R {input.fasta} \
        --intervals {input.bed} \
        -I {input.bam} \
        -BQSR {input.table} \
        -o {output}
        """

# Create table for recalibration

rule table_for_base_recal:
    input:
        bam="output/bam_recal/{sample}_dupmarked_realigned.bam",
        fasta=fasta_path,
        sites=known_sites_path,
        bed=bed_file_path
    output:
        temp("output/bam_recal/{sample}_dupmarked_realigned_recal.table")
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T BaseRecalibrator \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -knownSites {input.sites} \
        -o {output}
        """

# Realign BAM files by indels

rule indel_realigner:
    input:
        bam="output/bam_recal/{sample}_dupmarked_reheader.bam",
        fasta=fasta_path,
        indels=known_indels_path,
        bed=bed_file_path,
        targets="output/bam_recal/{sample}_dupmarked_realigner.intervals"
    output:
       temp("output/bam_recal/{sample}_dupmarked_realigned.bam")
    params:
        gatk=gatk_path
    shell:
        """ 
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T IndelRealigner \
        --consensusDeterminationModel KNOWNS_ONLY \
        -LOD 0.4 \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -known {input.indels} \
        -targetIntervals {input.targets} \
        -o {output}
        """

# Create Realignment Targets

rule realigner_target_creator:
    threads: 5
    input:
        bam="output/bam_recal/{sample}_dupmarked_reheader.bam",
        fasta=fasta_path,
        indels=known_indels_path,
        bed=bed_file_path
    output:
        temp("output/bam_recal/{sample}_dupmarked_realigner.intervals")
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T RealignerTargetCreator \
        -nt {threads} \
        -R {input.fasta} \
        -I {input.bam} \
        --intervals {input.bed} \
        -known {input.indels} \
        -o {output}
        """

# Update headers and index bam files

rule update_header_and_index:
    input:
        "output/bam_recal/{sample}_dupmarked.bam"
    output:
        bam_output=temp("output/bam_recal/{sample}_dupmarked_reheader.bam"),
        bam_index=temp("output/bam_recal/{sample}_dupmarked_reheader.bam.bai")
    params:
        header = lambda w: f"'s,^@RG.*,@RG\\tID:{w.sample}\\tSM:{w.sample}\\tLB:None\\tPL:Illumina,g'"
    shell:
        """
        samtools view -H {input} | \
        sed {params.header} | \
        samtools reheader - {input} > {output.bam_output}

        samtools index {output.bam_output} 
        """

# Mark duplicates

rule mark_duplicates:
    input: 
        bam="output/sorted/{sample}_sorted.bam",
        bam_index="output/sorted/{sample}_sorted.bam.bai"
    output: 
        dup_marked=temp("output/bam_recal/{sample}_dupmarked.bam"),
        metrics_file=temp("output/bam_recal/{sample}_picard_metrics_file.txt")
    params:
        picard=picard_path
    shell:
        """
        java -Djava.iodir=1000m -jar {params.picard} \
        MarkDuplicates AS=TRUE VALIDATION_STRINGENCY=LENIENT \
        I={input.bam} \
        O={output.dup_marked} \
        M={output.metrics_file}
        """
        
# Sort BAM for recalibration & index

rule index_sorted_bam: 
    input: 
        "output/sorted/{sample}_sorted.bam"
    output:
        temp("output/sorted/{sample}_sorted.bam.bai") 
    shell:
        "samtools index {input}"

rule samtools_sort:
    threads: 5
    input: 
        "output/mapped_reads/{sample}.bam"
    output:
        temp("output/sorted/{sample}_sorted.bam")
    shell:
        "samtools sort -@ {threads} {input} > {output}"

# Map reads 

rule bwa_map:
    threads: 5
    input:
        fasta_path,
        f"{source_dir}/{{sample}}_1.fastq.gz",
        f"{source_dir}/{{sample}}_2.fastq.gz"
    output:
        "output/mapped_reads/{sample}.bam"
    params:
        header = lambda w: f"@RG\\\\tID:{w.sample}\\\\tPL:ILLUMINA"
    shell:
         "bwa mem -t {threads} -M -R {params.header} {input} | samtools view -u -S - | samtools sort -n -o {output}"
