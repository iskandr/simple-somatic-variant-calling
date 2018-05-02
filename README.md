# simple-somatic-variant-calling
A simple genomics "pipeline" implemented via a single shell script. Aligns tumor & normal DNA sequencing data, marks duplicate reads in BAM files, and runs Strelka2 as a somatic variant caller. 

Expected input FASTQ files for normal and tumor DNA sequencing data. The pipeline will generate (in whichever directory it's invoked) a collection of BAM files and ultimately two VCFs: somatic SNVs and somatic indels. 

## Requirements

 * [sambamba](http://lomereiter.github.io/sambamba/)
 * [bamtools](https://github.com/pezmaster31/bamtools)
 * [Picard](https://broadinstitute.github.io/picard/)
 * [BWA](http://bio-bwa.sourceforge.net/)
 * [Strelka2](https://github.com/Illumina/strelka/)

## Invocation 
```sh
./call-somatic-variants $PATH_TO_NORMAL $NORMAL_FASTQ_PREFIX $PATH_TO_TUMOR $TUMOR_FASTQ_PREFIX
```

Generates BAMs and VCF file in the directory it was invoked. You can edit the shell script to specify a directory with an indexed human reference, otherwise the reference gets downloaded and indexed in the same directory. 
