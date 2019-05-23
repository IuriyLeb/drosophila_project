REFERENCE = "reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa"
WILDSTRAINS = ["GER1800140_Cs_LP180329001", "GER1800136_Ber_LP180329001", "GER1800138_OR_LP180329001"]
MUTANTS = ["GER1800139_X1_LP180329001", "GER1800137_ts3_LP180329001"]

rule all:
    input:
        expand("snpsift/{mutants}/{impact}_impact.vcf", mutants=MUTANTS, impact=['high', 'low', 'moderate', 'modifier']),
        expand("qc/fastqc_rep/{mutants}/report_{mutants}.html", mutants=MUTANTS),
        expand("qc/fastqc_rep/{mutants}/report_{mutants}.zip",  mutants=MUTANTS),
        expand("qc/fastqc_rep/{wildstrains}/report_{wildstrains}.html", wildstrains=WILDSTRAINS),
        expand("qc/fastqc_rep/{wildstrains}/report_{wildstrains}.zip", wildstrains=WILDSTRAINS),
        expand("snpeff/{mutants}/{mutants}_filtered_indels.vcf", mutants=MUTANTS),
        expand("snpeff/{mutants}/{mutants}_filtered_snp.vcf", mutants=MUTANTS)

rule snpsift_modifier:
        input:
            "snpeff/{mutants}/{mutants}_ann_chrx.vcf",
        output:
            "snpsift/{mutants}/modifier_impact.vcf"
        shell:
            '''java -Xmx4g -jar ~/snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'MODIFIER'" {input} > {output}'''

rule snpsift_low:
    input:
        "snpeff/{mutants}/{mutants}_ann_chrx.vcf",
    output:
        "snpsift/{mutants}/low_impact.vcf"
    shell:
        '''java -Xmx4g -jar ~/snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" {input} > {output}'''

rule snpsift_moderate:
    input:
        "snpeff/{mutants}/{mutants}_ann_chrx.vcf",
    output:
        "snpsift/{mutants}/moderate_impact.vcf"
    shell:
        '''java -Xmx4g -jar ~/snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" {input} > {output}'''

rule snpsift_high:
    input:
        "snpeff/{mutants}/{mutants}_ann_chrx.vcf",
    output:
        "snpsift/{mutants}/high_impact.vcf"
    shell:
        '''java -Xmx4g -jar ~/snpEff/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" {input} > {output}'''


rule snpeff:
    input:
      "results/{mutants}_chrx.vcf",
    output:
      "snpeff/{mutants}/{mutants}_ann_chrx.vcf",
    shell:
      '''java -Xmx4g -jar ~/snpEff/snpEff.jar -v BDGP6.86 {input} > {output};
      mv snpEff_genes.txt snpeff/{wildcards.mutants}/snpEff_genes.txt;
      mv snpEff_summary.html snpeff/{wildcards.mutants}/snpEff_summary.html'''

rule cut_x_chrom:
    input:
        vars = "snpeff/{mutants}/{mutants}_ann.vcf",
    output:
        x_chrom = "results/{mutants}_chrx.vcf"
    shell:
        "bcftools view {input.vars} --regions X -o {output.x_chrom}"

rule filter_snpsift_indels:
    input:
        "snpeff/{mutants}/{mutants}_indels.vcf"
    output:
        "snpeff/{mutants}/{mutants}_filtered_indels.vcf"
    shell:
        '''java -jar SnpSift.jar filter "((ReadPosRankSum > -20.0) & (QD > 2.0) & (FS < 200.0) & (SOR < 10.0))" {input} > {output}'''


rule filter_snpsift_snp:
    input:
        "snpeff/{mutants}/{mutants}_snp.vcf"
    output:
        "snpeff/{mutants}/{mutants}_filtered_snp.vcf"
    shell:
        '''java -jar SnpSift.jar filter "((ReadPosRankSum > -8.0) & (MQRankSum > -2.5) & (SOR < 3.0) & (QD > 2.0) & (FS < 60.0) & (MQ > 50.0))" {input} > {output}'''


rule filter_indels:
    input:
        "snpeff/{mutants}/{mutants}_ann.vcf"
    output:
        "snpeff/{mutants}/{mutants}_indels.vcf"
    shell:
        "vcftools —vcf {input} —keep-only-indels —out {output} —recode —recode-INFO-all"

rule filter_snp:
    input:
        "snpeff/{mutants}/{mutants}_ann.vcf"
    output:
        "snpeff/{mutants}/{mutants}_snp.vcf"
    shell:
        "vcftools —vcf {input} —remove-indels —out {output} —recode —recode-INFO-all"

rule snpeff_mut:
    input:
        "vcfeval/{mutants}/fp.vcf"
    output:
        "snpeff/{mutants}/{mutants}_ann.vcf",
    shell:
        '''java -Xmx4g -jar ~/snpEff/snpEff.jar -v BDGP6.86 {input} > {output};
        mv snpEff_genes.txt snpeff/{wildcards.mutants}/snpEff_genes.txt;
        mv snpEff_summary.html snpeff/{wildcards.mutants}/snpEff_summary.html'''

rule unzip_vars:
    input:
        "vcfeval/{mutants}/fp.vcf.gz"
    output:
        "vcfeval/{mutants}/fp.vcf"
    shell:
        "sudo gunzip -k {input}"


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
        sam="mutants/sam/{mutants}.sam",
        idx="reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fai"
    output:
        "mutants/bam/{mutants}.bam",
    threads: 4
    shell:
        "samtools view -@ {threads} -S -b {input.sam} > {output}"

rule index_for_gatk:
    input:
        "reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa"
    output:
        "reference/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fai"
    shell:
        "samtools faidx {input}"

rule bwa_wildstrains:
    input:
        REFERENCE,
        "wildstrains/trimmed/{wildstrains}_R1.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.paired.fastq.gz",
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

rule trimmommatic_wildstrains:
    input:
        'wildstrains/{wildstrains}_R1.fastq.gz',
        'wildstrains/{wildstrains}_R2.fastq.gz',
    output:
        "wildstrains/trimmed/{wildstrains}_R1.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R1.unpaired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.paired.fastq.gz",
        "wildstrains/trimmed/{wildstrains}_R2.unpaired.fastq.gz",
    threads: 4
    shell:
        "TrimmomaticPE -threads {threads} {input} {output} LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:30"


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
    wrapper:
        "0.34.0/bio/fastqc"

rule fastqc_wildstrains:
    input:
        "wildstrains/{wildstrains}_R1.fastq.gz"
    output:
        html="qc/fastqc_rep/{wildstrains}/report_{wildstrains}.html",
        zip="qc/fastqc_rep/{wildstrains}/report_{wildstrains}.zip"
    wrapper:
        "0.34.0/bio/fastqc"
