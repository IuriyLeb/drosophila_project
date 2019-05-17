REFERENCE = "reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa"
WILDSTRAINS = ["GER1800140_Cs_LP180329001", "GER1800136_Ber_LP180329001", "GER1800138_OR_LP180329001"]
MUTANTS = ["GER1800139_X1_LP180329001", "GER1800137_ts3_LP180329001"]

rule all:
    input:
        expand("snpeff/{mutants}/{mutants}_ann.vcf.gz", mutants=MUTANTS),
        # expand("qc/fastqc_rep/{mutants}/report_{mutants}.html", mutants=MUTANTS),
        # expand("qc/fastqc_rep/{mutants}/report_{mutants}.zip",  mutants=MUTANTS),
        # expand("qc/fastqc_rep/{wildstrains}/report_{wildstrains}.html", wildstrains=WILDSTRAINS),
        # expand("qc/fastqc_rep/{wildstrains}/report_{wildstrains}.zip", wildstrains=WILDSTRAINS)

rule snpeff:
    input:
      "results/{mutants}_chrx.vcf",
    output:
      "snpeff/{mutants}/{mutants}_ann.vcf.gz",
    shell:
      "java -Xmx4g -jar /home/iuriy/snpEff/snpEff.jar -v BDGP6.86 {input} > {output}"

rule cut_x_chrom:
    input:
        bed = "reference/xchrom.bed",
        vars = "vcfeval/{mutants}/fp.vcf",
    output:
        x_chrom = "results/{mutants}_chrx.vcf"
    shell:
        "bedtools intersect -b {input.bed} -a {input.vars} > {output.x_chrom}"

rule unzip_vars:
    input:
        "vcfeval/{mutants}/fp.vcf.gz"
    output:
        "vcfeval/{mutants}/fp.vcf"
    shell:
        "gunzip {input}"


rule vcfeval:
    input:
        ref = REFERENCE,
        wild="merged/wildstrains/wildstrains_merged.vcf.gz",
        mut="gatk/mutants/{mutants}.g.vcf.gz",
        dir_sdf="Drosophila_melanogaster.BDGP6.dna_sm.toplevel"
    params:
        outdir = "vcfeval/{mutants}/"
    output:
        "vcfeval/{mutants}/fp.vcf.gz",
    shell:
        "rm -r {params.outdir}; ~/rtg-tools/rtg vcfeval -b {input.mut} -c {input.wild} -o {params.outdir} -t {input.dir_sdf}"

rule make_sdf:
    input:
        ref = REFERENCE
    params:
        outdir= "Drosophila_melanogaster.BDGP6.dna_sm.toplevel",
    output:
        "Drosophila_melanogaster.BDGP6.dna_sm.toplevel/done"
    shell:
        "~/rtg-tools/rtg format -o {params.outdir} {input.ref}"

rule mergewt:
    input:
      expand("gatk/wildstrains/{wildstrains}.g.vcf.gz", wildstrains = WILDSTRAINS)
    output:
      "merged/wildstrains/wildstrains_merged.vcf.gz"
    shell:
      "gatk MergeVcfs -I {input} -O {output}"

rule gatk_mutants:
    input:
        "picard/mutants/{mutants}_picard.bam",
    output:
        "gatk/mutants/{mutants}.g.vcf.gz"
    threads: 4
    shell:
        "gatk --java-options -Xmx3g HaplotypeCaller -R {REFERENCE} -I {input} -O {output}"

rule gatk_wildstrains:
    input:
        "picard/wildstrains/{wildstrains}_picard.bam",
    output:
        "gatk/wildstrains/{wildstrains}.g.vcf.gz",
    threads: 4
    shell:
        "gatk --java-options -Xmx3g HaplotypeCaller -R {REFERENCE} -I {input} -O {output}"

rule index_picard_wildstrains:
    input:
        "picard/wildstrains/{wildstrains}_picard.bam"
    output:
        "picard/wildstrains/{wildstrains}_picard.bam.bai"
    threads: 4
    shell:
        "samtools index {input} {output}"

rule index_picard_mutants:
    input:
        "picard/mutants/{mutants}_picard.bam"
    output:
        "picard/mutants/{mutants}_picard.bam.bai"
    threads: 4
    shell:
        "samtools index {input} {output}"

rule picard_wildstrains:
    input:
        "wildstrains/bam_s/{wildstrains}_sorted.bam"
    output:
        "picard/wildstrains/{wildstrains}_picard.bam"
    threads: 4
    shell:
        "PicardCommandLine AddOrReplaceReadGroups I={input} O={output} RGID=1 RGLB=libraryname RGPL=illumina RGPU=unitname RGSM=20"

rule picard_mutants:
    input:
        "mutants/bam_s/{mutants}_sorted.bam"
    output:
        "picard/mutants/{mutants}_picard.bam"
    threads: 4
    shell:
        "PicardCommandLine AddOrReplaceReadGroups I={input} O={output} RGID=1 RGLB=libraryname RGPL=illumina RGPU=unitname RGSM=20"

rule samtools_ind_wildstrains:
    input:
        "wildstrains/bam_s/{wildstrains}_sorted.bam"
    output:
        "wildstrains/bam_s/{wildstrains}_sorted.bam.bai"
    threads: 4
    shell:
        "samtools index {input} {output}"

rule samtools_ind_mutants:
    input:
        "mutants/bam_s/{mutants}_sorted.bam"
    output:
        "mutants/bam_s/{mutants}_sorted.bam.bai"
    threads: 4
    shell:
        "samtools index {input} {output}"

rule samtools_sort_wildstrains:
    input:
        "wildstrains/bam/{wildstrains}.bam"
    output:
        "wildstrains/bam_s/{wildstrains}_sorted.bam"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools_sort_mutants:
    input:
        "mutants/bam/{mutants}.bam",
    output:
        "mutants/bam_s/{mutants}_sorted.bam"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools_wildstrains:
    input:
        "wildstrains/sam/{wildstrains}.sam",
    output:
        "wildstrains/bam/{wildstrains}.bam",
    threads: 4
    shell:
        "samtools view -@ {threads} -S -b {input} > {output}"

rule samtools_mutants:
    input:
        "mutants/sam/{mutants}.sam",
    output:
        "mutants/bam/{mutants}.bam",
    threads: 4
    shell:
        "samtools view -@ {threads} -S -b {input} > {output}"

rule bwa_wildstrains:
    input:
        REFERENCE,
        "wildstrains/trimmed/{wildstrains}.1.fastq.gz",
        "wildstrains/trimmed/{wildstrains}.2.fastq.gz",
    output:
        "wildstrains/sam/{wildstrains}.sam",
    threads: 4
    shell:
        "bwa mem -t {threads} {input} > {output}"

rule bwa_mutants:
    input:
        REFERENCE,
        "mutants/trimmed/{mutants}_R1.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.paired.fastq.gz",
    output:
        "mutants/sam/{mutants}.sam",
    threads: 4
    shell:
        "bwa mem -t {threads} {input} > {output}"

# rule trimmommatic_wildstrains:
#     input:
#         'wildstrains/{wildstrains}_R1.fastq.gz',
#         'wildstrains/{wildstrains}_R2.fastq.gz',
#     output:
#         "wildstrains/trimmed/{wildstrains}_R1.paired.fastq.gz",
#         "wildstrains/trimmed/{wildstrains}_R1.unpaired.fastq.gz",
#         "wildstrains/trimmed/{wildstrains}_R2.paired.fastq.gz",
#         "wildstrains/trimmed/{wildstrains}_R2.unpaired.fastq.gz",
#     threads: 4
#     shell:
#         "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"

rule trimmomatic_pe_wild:
    input:
        r1="wildstrains/{wildstrains}_R1.fastq.gz",
        r2="wildstrains/{wildstrains}_R2.fastq.gz",
    output:
        r1="wildstrains/trimmed/{wildstrains}.1.fastq.gz",
        r2="wildstrains/trimmed/{wildstrains}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="wildstrains/trimmed/{wildstrains}.1.unpaired.fastq.gz",
        r2_unpaired="wildstrains/trimmed/{wildstrains}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{wildstrains}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["LEADING:20", "TRAILING:20", "SLIDINGWINDOW:6:20", "MINLEN:30"],
    threads:
        4
    wrapper:
        "0.34.0/bio/trimmomatic/pe"


rule trimmommatic_mutants:
    input:
        "mutants/{mutants}_R1.fastq.gz",
        "mutants/{mutants}_R2.fastq.gz",
    output:
        "mutants/trimmed/{mutants}_R1.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R1.unpaired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.unpaired.fastq.gz",
    threads: 4
    shell:
        "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"


rule fastqc_mutants:
    input:
        "mutants/{mutants}_R1.fastq.gz"
    output:
        html="qc/fastqc_rep/{mutants}/report_{mutants}.html",
        zip="qc/fastqc_rep/{mutants}/report_{mutants}.zip"
    log:
        "logs/fastqc/{mutants}.log"
    wrapper:
        "0.34.0/bio/fastqc"

rule fastqc_wildstrains:
    input:
        "wildstrains/{wildstrains}_R1.fastq.gz"
    output:
        html="qc/fastqc_rep/{wildstrains}/report_{wildstrains}.html",
        zip="qc/fastqc_rep/{wildstrains}/report_{wildstrains}.zip"
    log:
        "logs/fastqc/{wildstrains}.log"
    wrapper:
        "0.34.0/bio/fastqc"
