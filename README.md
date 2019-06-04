# Whole Genome Drosophila Sequence
## Description
Drosophila(*D.melanogaster*) is a widely known and popular model organism. Due to the small size of the Drosophila genome and its rapid reproduction, it is fairly easy to reproduce various conditions and diseases with their subsequent study. In this work, two Drosophila lines (X1 and ts3) with disturbances in the structure and functioning of the nervous system, which were compared with three wild-type lines (Canton-S, Berlin, Oregon-R), without physiological pathologies, were studied. Both of these lines can potentially be used to study Williams syndrome in humans.
## Goals and objectives
Our goal was to detect mutations responsible for the formation of mutant phenotypes.
To achieve this goal the following tasks were set.
* Basic reads quality control
* Find variants for all lines
* Find out the differences between the genomes of wild and mutant lines
* Identify mutations responsible for the formation of mutant phenotypes.  
**All stages are included in the general pipeline as a Snakefile file, presented in this repository.**
**Some rules require `sudo` privileges, be ready to enter your super-user password**
## Methods
### Necessary software installing.
* FastQC(v0.11.5) - `sudo apt install fastqc`
* TrimmomaticPE(v0.22) - `sudo apt install trimmomatic`
* BWA(v0.7.17-r1188) - `sudo apt install bwa`
* Samtools(v1.7) - `sudo apt install samtools`
* Picard(v1.138) - `sudo apt install picard`
* [GATK](https://software.broadinstitute.org/gatk/download/index)(v4.1.0.0)
* [RTGTools](https://github.com/RealTimeGenomics/rtg-tools)(v3.10.1) - must be available in the following path: `~/rtg-tools/rtg`
* Bedtools(v2.26.0) - `sudo apt install bedtools`
* vcftools(v0.1.17) - `sudo apt install vcftools`
* [SNPEff](http://snpeff.sourceforge.net/download.html)(v4.3t) - must be available in the following path: `~/snpEff/snpEff.jar`
  * Final SnpSift filtration runs with following parameters:
    * `"((ReadPosRankSum > -8.0) & (MQRankSum > -2.5) & (SOR < 3.0) & (QD > 2.0) & (FS < 60.0) & (MQ > 50.0))"` for SNP's
    * `"((ReadPosRankSum > -20.0) & (QD > 2.0) & (FS < 200.0) & (SOR < 10.0))"` for Indels
## Files description    
1. Snakefile - executable file for Snakemake. To run, simply type `snakemake` in terminal in project directory.   
2. `pipeline_*.pdf` - pipeline scheme.
  * `pipeline_wo_fqc.pdf` - pipeline scheme without fastqc rule added(more compact).
  * `pipeline_with_fqc.pdf` - pipeline scheme with fastqc rule added(big picture).
3. `snpEff` - directory with snpEff results and summary files.
4. `/files` - directory with result images.

## Results
1. We obtained variants for each of the mutant and wild-type lines. Then we summed up wild-type variants and then duducted them from the mutant lines (X1, ts3) variants with vcfeval:   

|                | ts3 | x1 |
|----------------|-----|----|
| before vcfeval |   815608  |  803690  |
| after vcfeval  |  326780   |  302019  |
| in chrX        |   46165  |   46032 |

![](https://github.com/IuriyLeb/drosophila_project/blob/master/files/variants.jpg?raw=true)
### Common variants for ts3 and x1
<img src="https://github.com/IuriyLeb/drosophila_project/raw/master/files/common.jpg?raw=true" alt="" width="450" height="300">
2. These variations are annotated and filtered according to their effect.

![](https://github.com/IuriyLeb/drosophila_project/blob/master/files/snp_indels.jpg?raw=true)

3. The filtered variations are divided into groups by significance for manual analysis.
4. The entire pipeline is made in Snakemake
