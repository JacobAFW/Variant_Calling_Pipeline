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

# Create bed for each chromosome/contig - to be used for sub-setting jobs (eg joint calling)

with open(config['bed_file_path']) as f:
    chrom_bed = f.read().splitlines()
    chrom_bed = [p.split('\t') for p in chrom_bed]

for i in chrom_bed:
    content = '\t'.join(i)
    file_name = i[0]
    file_path = 'data/chromosomes/' + file_name + '.bed'
    f = open(file_path, 'w') # amend this line - how to amend bed file path and use chrom as file name
    f.write(content)

f.close()

# Define final files

rule all:
    input:
        f"{fasta_path}.fai",
        f"{fasta_path}.dict",
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
        expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample = SAMPLES),
        "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz",
        expand("output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz", chromosome = CHROMOSOME), 
        expand("output/calling/bcftools/input_bam_files.list"),
        expand("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz", chromosome = CHROMOSOME),
        expand("output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz.tbi", chromosome = CHROMOSOME),
        expand("output/calling/consensus/{chromosome}_consensus.vcf.gz", chromosome = CHROMOSOME),
        expand("output/calling/consensus/{chromosome}_consensus.vcf.gz.tbi", chromosome = CHROMOSOME),
        expand("output/calling/consensus/{chromosome}.txt", chromosome = CHROMOSOME),
        "output/calling/consensus/Consensus.vcf.gz",
        "output/calling/consensus/Consensus.vcf.gz.tbi"

# Define local rules - not run with scheduler
localrules: 
    all, bam_input_list, index_ref

# Concatenate chromosome-based consensus VCFs 
rule concat_vcfs:
    input:
        expand("output/calling/consensus/{chromosome}_consensus.vcf.gz", chromosome = CHROMOSOME)
    params:
        vcf = lambda w: " ".join(expand("output/calling/consensus/{chromosome}_consensus.vcf.gz", chromosome = CHROMOSOME))
    output:
        vcf="output/calling/consensus/Consensus.vcf.gz",
        tbi="output/calling/consensus/Consensus.vcf.gz.tbi"
    shell:
        """
        bcftools concat -o {output.vcf} {params.vcf}
        bcftools index -t -o {output.tbi} {output.vcf}
        """

# Take a consenus of GATK and bcftools
rule consensus_of_vcfs:
    input:
        bcftools="output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz",
        gatk="output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz"
    output:
        txt="output/calling/consensus/{chromosome}.txt",
        vcf="output/calling/consensus/{chromosome}_consensus.vcf.gz",
        tbi="output/calling/consensus/{chromosome}_consensus.vcf.gz.tbi"
    params:
        header = lambda w: f"'%CHROM\\t%POS\\n'"
    shell:
        """
        bcftools query -f {params.header} {input.bcftools} > {output.txt}
        bcftools filter -R {output.txt} -o {output.vcf} {input.gatk}
        bcftools index -t -o {output.tbi} {output.vcf}
        """

# Run bcftools for each chromosome
rule bcftools_caller:
    input:
        input_bam_files="output/calling/bcftools/input_bam_files.list",
        fasta=fasta_path,
        bed="data/chromosomes/{chromosome}.bed"
    output:
        vcf="output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz",
        tbi="output/calling/bcftools/bcftools_genotyped_{chromosome}.vcf.gz.tbi"
    shell:
        """
        bcftools mpileup --threads 2 -f {input.fasta} -b {input.input_bam_files} -R {input.bed} | bcftools call --threads 2 -m -Oz -a FORMAT/GQ,FORMAT/GP,INFO/PV4 -v -o {output.vcf}
        bcftools index --threads 2 -t -o {output.tbi} {output.vcf}
        """

# Create input list of bam files for bcftools

rule bam_input_list:
    input:
        bam=expand("output/bam_recal/{sample}_recalibrated.bam", sample = SAMPLES)
    output:
        temp("output/calling/bcftools/input_bam_files.list")
    run:
        import glob

        bam_list = glob.glob('output/bam_recal/*_recalibrated.bam')
        #bam_list = [sub.replace('output/bam_recal/', '') for sub in bam_list] 

        file = open('output/calling/bcftools/input_bam_files.list', 'w')
        for item in bam_list:
            file.write(item+"\n")

        file.close()

# Joint-call variants

rule joint_genotyping:
    input:
        vcf="output/calling/gatk/gvcf/GATK_combined.g.vcf.gz",
        fasta=fasta_path,
        bed="data/chromosomes/{chromosome}.bed"
    output:
        "output/calling/gatk/joint/gatk_genotyped_{chromosome}.vcf.gz"
    params:
        gatk=gatk_path
    shell:
        """
        java -Djava.iodir=1000m -Xms3200m -Xmx3600m -jar {params.gatk} \
        -T GenotypeGVCFs \
        -nt 3 \
        -R {input.fasta} \
        -L {input.bed} \
        -V {input.vcf} \
        -o {output}
        """
            
# Combine gVCFs

rule combine_gvcfs:
    input:
        fasta=fasta_path,
        gvcf=expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample = SAMPLES)
    output:
        "output/calling/gatk/gvcf/GATK_combined.g.vcf.gz"
    params:
        gatk=gatk_path,
        gvcfs = lambda w: " -V " + " -V ".join(expand("output/calling/gatk/gvcf/{sample}.g.vcf.gz", sample = SAMPLES))
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
        "output/calling/gatk/gvcf/{sample}.g.vcf.gz"
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

rule samtools_sort:
    threads: 5
    input: 
        "output/mapped_reads/{sample}.bam"
    output:
        bam_output=temp("output/sorted/{sample}_sorted.bam"),
        bam_index=temp("output/sorted/{sample}_sorted.bam.bai") 
    shell:
        """
        samtools sort -@ {threads} {input} > {output.bam_output}

        samtools index {output.bam_output}
        """

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

# Index reference genome

rule index_ref:
    input:
        fasta_path
    output:
        index=f"{fasta_path}.fai",
        dictionary=f"{fasta_path}.dict"
    params:
        picard_path
    shell:
        """
        samtools faidx -o {output.index} {input}
        java -jar {params} CreateSequenceDictionary R={input} O={output.dictionary}
        """