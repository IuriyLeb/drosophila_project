REFERENCE = "reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa"
WILDSTRAINS = ["GER1800140_Cs_LP180329001", "GER1800136_Ber_LP180329001", "GER1800138_OR_LP180329001"]
MUTANTS = ["GER1800139_X1_LP180329001", "GER1800137_ts3_LP180329001"]
MERGED = ["GER1800140_Cs_LP180329001", "GER1800136_Ber_LP180329001", "GER1800138_OR_LP180329001", "GER1800139_X1_LP180329001", "GER1800137_ts3_LP180329001"]

rule all:
    input:
        expand("vcfeval/vcfeval_ann_{mutants}.vcf.gz", mutants=MUTANTS)

rule snpeff:
    input:
      "vcfeval/{mutants}/fp.vcf.gz",
    output:
      "vcfeval/{mutants}/vcfeval_ann_{mutants}.vcf.gz",
    shell:
      "java -Xmx4g -jar /home/iuriy/snpEff/snpEff.jar -v BDGP6.86 {input} > {output}"


rule vcfeval:
    input:
        ref = REFERENCE,
        wild="merged/wildstrains/wildstrains_merged.vcf.gz",
        mut="gatk/mutants/{mutants}.g.vcf.gz",
        inp_dir="Drosophila_melanogaster.BDGP6.dna_sm.toplevel/"
    output:
        outdir="vcfeval/{mutants}/",
    shell:
        "~/rtg-tools/rtg vcfeval -b {input.mut} -c {input.wild} -o {{output.outdir}} -t {input.inp_dir}

rule make_sdf:
    input:
        ref = REFERENCE
    output:
        outdir="Drosophila_melanogaster.BDGP6.dna_sm.toplevel/"
    shell:
        "~/rtg-tools/rtg format -o {output.outdir} {input.ref}"

rule mergewt:
    input:
      cs = 'gatk/wildstrains/{{}}.g.vcf.gz'.format(WILDSTRAINS[0]),
      ber = 'gatk/wildstrains/{{}}.g.vcf.gz'.format(WILDSTRAINS[1]),
      ore = 'gatk/wildstrains/{{}}.g.vcf.gz'.format(WILDSTRAINS[2]),
    output:
      "merged/wildstrains/wildstrains_merged.vcf.gz"
    shell:
      "gatk MergeVcfs -I {input.cs} -I {input.ber} -I {input.ore} -O {output}"

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
        "picard//mutants/{mutants}_picard.bam.bai"
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
        "mutants/bam/{mutants}.bam"
    output:
        "mutants/bam_s/{mutants}_sorted.bam"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools_wildstrains:
    input:
        "wildstrains/sam/{wildstrains}.sam"
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
        "wildstrains/trimmed/{wildstrains}_R1.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.paired.fastq.gz"
    output:
        "wildstrains/sam/{wildstrains}.sam",
    threads: 4
    shell:
        "bwa mem -t {threads} {input} > {output}"

rule bwa_mutants:
    input:
        REFERENCE,
        "mutants/trimmed/{mutants}_R1.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.paired.fastq.gz"
    output:
        "mutants/sam/{mutants}.sam",
    threads: 4
    shell:
        "bwa mem -t {threads} {input} > {output}"

rule trimmommatic_wildstrains:
    input:
        "wildstrains/{wildstrains}_R1.fastq.gz",
        "wildstrains/{wildstrains}_R2.fastq.gz"
    output:
        "wildstrains/trimmed/{wildstrains}_R1.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R1.unpaired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.unpaired.fastq.gz"
    threads: 4
    shell:
        "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"

rule trimmommatic_mutants:
    input:
        "mutants/{mutants}_R1.fastq.gz",
        "mutants/{mutants}_R2.fastq.gz"
    output:
        "mutants/trimmed/{mutants}_R1.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R1.unpaired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.paired.fastq.gz",
        "mutants/trimmed/{mutants}_R2.unpaired.fastq.gz"
    threads: 4
    shell:
        "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"
