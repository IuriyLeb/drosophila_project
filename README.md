# Whole Genome Drosophila Sequence
## Description
Drosophila(*D.melanogaster*) is a widely known and popular model organism. Due to the small size of the Drosophila genome and its rapid reproduction, it is fairly easy to reproduce various conditions and diseases with their subsequent study. In this work, two Drosophila lines (X1 and ts3) with disturbances in the structure and functioning of the nervous system, which were compared with three wild-type lines, without physiological pathologies, were studied. Both of these lines can potentially be used to study Williams syndrome in humans.
## Goals and objectives
Our goal was to detect mutations responsible for the formation of mutant phenotypes.
To achieve this goal the following tasks were set.
* Basic reads quality control
* Find variants for all lines
* Find out the differences between the genomes of wild and mutant lines
* Identify mutations responsible for the formation of mutant phenotypes.  
**All stages are included in the general pipeline as a Snakefile file, presented in this repository.**
**Some rules require `sudo` privileges**
## Methods
### Necessary software installing.
* FastQC - `sudo apt install fastqc`
* TrimmomaticPE - `sudo apt install trimmomatic`
* BWA - `sudo apt install bwa`
* Samtools - `sudo apt install samtools`
* Picard - `sudo apt install picard`
* [GATK](https://software.broadinstitute.org/gatk/download/index)
* [RTGTools](https://github.com/RealTimeGenomics/rtg-tools) - must be available in the following path: `~/rtg-tools/rtg`
* Bedtools - `sudo apt install bedtools`
* [SNPEff](http://snpeff.sourceforge.net/download.html) - must be available in the following path: `~/snpEff/snpEff.jar`
## Files description
1. Snakefile - executable file for Snakemake. To run, simply enter `snakemake` in project directory
2. `./reference/xchrom.bed` - additional file for intersect with `.vcf` to leave only one chromosome.
3. `dag.pdf` - pipeline scheme
## Results
