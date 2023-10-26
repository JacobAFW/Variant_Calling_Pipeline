# Snakefile 

# Prerequisites (on Gadi): source /g/data/pq84/malaria/snakemake_pipeline/snakemake_setup.sh

# Define input files - defining the samples to be used for the {sample} wildcard

SAMPLES1, = glob_wildcards("data/{sample}_1.fastq.gz")
SAMPLES2, = glob_wildcards("data/{sample}_2.fastq.gz")


# Define final files

rule all:
    input:
        expand("mapped_reads/{sample}.bam", sample = SAMPLES1)
        expand("mapped_reads/{sample}.bam.bai", sample = SAMPLES1)
        expand("sorted/{sample}_sorted.bam", sample = SAMPLES1)
        expand("sorted/{sample}_sorted.bam.bai", sample = SAMPLES1)
        expand("bam_recal/{sample}_dupmarked.bam", sample = SAMPLES1)
        expand("bam_recal/{sample}_picard_metrics_file.txt", sample = SAMPLES1)
        expand("bam_recal/{sample}_dupmarked_reheader.bam", sample = SAMPLES1)
        expand("bam_recal/{sample}_dupmarked_realigner_intervals", sample = SAMPLES1)
        expand("bam_recal/{sample}_dupmarked_realigned_bam", sample = SAMPLES1)


# Realign by BAM files by indels
rule indel_realigner:
    input:
        bam="bam_recal/{sample}_dupmarked_reheader.bam",
        ref_genome_fasta="/data/Knowlesijwestaway/ref_genomes/PKA1H1/fasta/strain_A1_H.1.Icor.fasta"
        known_indels="/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/variants/TRUE_INDELs.vcf"
        bed_file="/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/strain_A1_H.1.Icor.fasta.bed"
        targets="bam_recal/{sample}_dupmarked_realigner_intervals"
    output:
        "bam_recal/{sample}_dupmarked_realigned_bam"
    shell:
        """ java -Djava.iodir=$PBS_JOBFS -Xms3200m -Xmx3600m -jar $GATK 
        -T IndelRealigner 
        --consensusDeterminationModel KNOWNS_ONLY 
        -LOD 0.4 
        -R $INDEXTDIR 
        -I {input.bam} 
        --intervals {input.bed_file}
        -known {input.known_indels} 
        -targetIntervals {input.targets}
        -o $OUTDIR/SAMPLE.dupmarked.realigned.bam
        """

# Create Realignment Targets
rule realigner_target_creator:
    input:
        bam="bam_recal/{sample}_dupmarked_reheader.bam",
        ref_genome_fasta="/data/Knowlesijwestaway/ref_genomes/PKA1H1/fasta/strain_A1_H.1.Icor.fasta"
        known_indels="/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/variants/TRUE_INDELs.vcf"
        bed_file="/g/data/pq84/malaria/Parasite_and_human_genetic_risk_factors_for_Pk_malaria/data/ref_genomes/PKA1H1/strain_A1_H.1.Icor.fasta.bed"
    output:
        "bam_recal/{sample}_dupmarked_realigner_intervals"
    shell:
        """java -Djava.iodir=$PBS_JOBFS -Xms3200m -Xmx3600m -jar $GATK 
        -T RealignerTargetCreator 
        -nt 5
        -R {input.ref_genome_fasta}
        -I {input.bam}
        --intervals {input.bed_file}
        -known {input.known_indels} 
        -o {output}
        """

# Bam Recalibration
rule update_header_and_index:
    input:
        "bam_recal/{sample}_dupmarked.bam"
    output:
        "bam_recal/{sample}_dupmarked_reheader.bam"
    shell:
    """samtools view -H {input} |
    sed 's,^@RG.*,@RG\tID:{sample}\tSM:{sample}\tLB:None\tPL:Illumina,g' |
    samtools reheader - {input} > {output}
    """

rule mark_duplicates:
    input: 
        "sorted/{sample}_sorted.bam",
        "sorted/{sample}_sorted.bam.bai"
    output: 
        temp(dup_marked="bam_recal/{sample}_dupmarked.bam"),
        temp(metrics_file="bam_recal/{sample}_picard_metrics_file.txt")
    shell:
        """java -Djava.iodir=$PBS_JOBFS -jar $PICARD
        MarkDuplicates AS=TRUE VALIDATION_STRINGENCY=LENIENT
        I={input}
        O={output.dup_marked}
        M={output.metrics_file}
        """
        
# Sort BAM for recalibration & index
rule index_sorted_bam: 
    input: 
        "sorted/{sample}_sorted.bam"
    output:
        temp("sorted/{sample}_sorted.bam.bai") 
    shell:
        "samtools index -o {output} {input}"

rule samtools_sort:
    input: 
        "mapped_reads/{sample}.bam"
    output:
        temp("sorted/{sample}_sorted.bam")
    shell:
        "samtools sort {input} > {output}"

# Map reads & index
rule index_bam:
    input: 
        "mapped_reads/{sample}.bam"
    output:
        "mapped_reads/{sample}.bam.bai"
    shell:
        "samtools index -o {output} {input}"

rule bwa_map:
    threads:5
    input:
        "/data/Knowlesijwestaway/ref_genomes/PKA1H1/fasta/strain_A1_H.1.Icor.fasta",
        "data/{sample}_1.fastq.gz",
        "data/{sample}_2.fastq.gz"
    output:
        "mapped_reads/{sample}.bam"
    shell:
         "bwa mem -t 5 -M -R '@RG\tID:{sample}\tPL:ILLUMINA' {input} | samtools view -u -S - | samtools sort -n -o {output}"

