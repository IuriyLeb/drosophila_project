
REFERENCE = "Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa"
WILDSTRAINS = ["GER1800140_Cs_LP180329001.g", "GER1800136_Ber_LP180329001.g", "GER1800138_OR_LP180329001.g"]

#Need to add bwa-index rule
rule snpEff:
    input:
      "wildstrains_merged.vcf.gz"
    output:
      "wildstrains_merged_ann.vcf.gz"
    shell:
      "java -Xmx4g -jar /home/iuriy/snpEff/snpEff.jar -v BDGP6.86 {input} > {output}"

rule mergewt:
    input:
      expand('{sample}.vcf.gz', sample=WILDSTRAINS)
    output:
      "wildstrains_merged.vcf.gz"
    shell:
      "vcf-merge {input} | bgzip > {output}"

rule gatk:
    input:
        ref="Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa",
        bam="{sample}_picard.bam",
        bai="{sample}_picard.bam.bai"
    output:
        "{sample}.g.vcf.gz"
    threads: 1
    shell:
        "gatk --java-options -Xmx3g HaplotypeCaller -R {input.ref} -I {input.bam} -O {output}"

rule picard:
    input:
        "{sample}_sorted.bam"
    output:
        "{sample}_picard.bam"
    threads: 1
    shell:
        "PicardCommandLine AddOrReplaceReadGroups I={input} O={output} RGID=1 RGLB=libraryname RGPL=illumina RGPU=unitname RGSM=20"

rule samtools_ind:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    threads: 1
    shell:
        "samtools index {input} {output}"

rule samtools_sort:
    input:
        "{sample}.bam"
    output:
        "{sample}_sorted.bam"
    threads: 1
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools:
    input:
        "{sample}.sam"
    output:
        "{sample}.bam",
    threads: 1
    shell:
        "samtools view -@ {threads} -S -b {input} > {output}"

rule bwa:
    input:
        REFERENCE,
        "{sample}_R1.paired.fastq.gz",
        "{sample}_R2.paired.fastq.gz"
    output:
        "{sample}.sam",
    threads: 1
    shell:
        "bwa mem -t {threads} {input} > {output}"

rule trimm:
    input:
        "{sample}_R1.fastq.gz",
        "{sample}_R2.fastq.gz"
    output:
        "{sample}_R1.paired.fastq.gz",
        "{sample}_R1.unpaired.fastq.gz",
        "{sample}_R2.paired.fastq.gz",
        "{sample}_R2.unpaired.fastq.gz"
    threads: 4
    shell:
        "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"


rule fastqc:
    input:
        "{sample}.fastq.gz",
    output:
        "{sample}_fastqc.zip",
        "{sample}_fastqc.html",
    threads: 1
    shell:
        "fastqc {input}"
