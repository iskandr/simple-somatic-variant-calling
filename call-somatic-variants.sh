#!/bin/bash
set -e

# machine configuration
NUMBER_PROCESSORS=32
MEMORY_LIMIT=32GB

# remote source for reference genome
REFERENCE_FASTA_SOURCE_SERVER=ftp://ftp.1000genomes.ebi.ac.uk
REFERENCE_FASTA_SOURCE_DIR=vol1/ftp/technical/reference/phase2_reference_assembly_sequence
REFERENCE_FASTA_SOURCE_NAME=hs37d5.fa.gz
REFERENCE_FASTA_SOURCE="$REFERENCE_FASTA_SOURCE_SERVER/$REFERENCE_FASTA_SOURCE_DIR/$REFERENCE_FASTA_SOURCE_NAME"

# local reference genome location
REFERENCE_DIR=.
REFERENCE_FASTA_NAME="hs37d5.fasta"
REFERENCE_FASTA_PATH="$REFERENCE_DIR/$REFERENCE_FASTA_NAME"
REFERENCE_INDEX_PATH="$REFERENCE_FASTA_PATH.bwt"

# normal sample
NORMAL_FASTQ_DIR=.
# first commandline argument is prefix for normal FASTQ files
NORMAL_FASTQ_PREFIX=${1:-'Z10-1002'}
NORMAL_SAMPLE_NAME="NORMAL"

# tumor sample
TUMOR_FASTQ_DIR=.
# second commandline argument is prefix for tumor FASTQ files
TUMOR_FASTQ_PREFIX=${2:-'Z10-1003'}
TUMOR_SAMPLE_NAME="TUMOR"

echo "Quick & Dirty Somatic Variant Calling Pipeline";
echo "---";
echo "Normal FASTQ location: $NORMAL_FASTQ_DIR/$NORMAL_FASTQ_PREFIX*.fastq";
echo "Tumor FASTQ location: $TUMOR_FASTQ_DIR/$TUMOR_FASTQ_PREFIX*.fastq";
echo "---";

function download_and_index_reference_genome() {
        echo "-- download_and_index_reference_genome";
        if [ ! -e $REFERENCE_FASTA_PATH ]; then
                echo "Couldn't find reference file $REFERENCE_FASTA_PATH";
                echo "Downloading from $REFERENCE_FASTA_SOURCE..."
                wget $REFERENCE_FASTA_SOURCE;
                echo "Decompressing downloaded reference genome..."
                time gunzip hs37d5.fa.gz;
                mv hs37d5.fa $REFERENCE_FASTA_PATH;
        else
                echo "Using reference: $REFERENCE_FASTA_PATH"
        fi;
        if [ ! -s $REFERENCE_INDEX_PATH ]; then
                echo "Creating index for $REFERENCE_FASTA_PATH"
                bwa index $REFERENCE_FASTA_PATH
        fi;
}

function align_fastq_pairs() {
        # 1) Align every FASTQ pair into multiple BAM files
        # 2) Merge them into single BAM file per sample
        local FASTQ_DIR=$1;
        local FASTQ_PREFIX=$2;
        local SAMPLE=$3;
        echo "-- align_fastq_pairs";
        echo "  FASTQ_DIR=$FASTQ_DIR";
        echo "  FASTQ_PREFIX=$FASTQ_PREFIX";
        echo "  SAMPLE=$SAMPLE";
        for R1_fastq in $FASTQ_DIR/$FASTQ_PREFIX*_R1*.gz ; do
                R2_fastq=`echo $R1_fastq | sed -e 's/_R1/_R2/g'`
                if [ ! -e $R2_fastq ]; then
                        echo "Couldn't find R2 ($R2_fastq) corresponding to $R1_fastq"
                        exit 1;
                fi;
                echo "R1: $R1_fastq";
                echo "R2: $R2_fastq";
                # make a local file name for the BAM we're going to generate from each FASTQ pair
                local READ_GROUP=`basename $R1_fastq | sed -e 's/_R1//g' | sed -e 's/.fastq.gz//g'`
                local BAM="$READ_GROUP.bam"
                # test if $BAM exists and is non-empty
                if [ ! -s $BAM ]; then
                        echo "Generating BAM file $BAM";
                        local BWA_COMMAND="bwa mem -M \
                                -t $NUMBER_PROCESSORS \
                                -R '@RG\tID:$READ_GROUP\tSM:$SAMPLE\tLB:$SAMPLE\tPL:ILLUMINA' \
                                $REFERENCE_FASTA_PATH \
                                $R1_fastq $R2_fastq | samtools view -Sb - > $BAM"
                        echo $BWA_COMMAND;
                        eval "time $BWA_COMMAND";
                        echo "---";
                else
                        echo "BAM file $BAM already exists";
                fi;
        done;
        echo "Merging all BAM files for $SAMPLE into single alignment file";
        bamtools merge -in $FASTQ_DIR/$FASTQ_PREFIX*.bam  -out $SAMPLE.bam;
}

function process_alignments() {
        # sort and index BAM file
        local $BAM=$1;
        local SORTED_BAM=`echo $BAM | sed -e 's/.bam/.sorted.bam'`;
        echo "Sorting $BAM to generate $SORTED_BAM";
        sambamba sort \
                --memory-limit $MEMORY_LIMIT \
                --show-progress \
                --nthreads $NUMBER_PROCESSORS \
                --out $SORTED_BAM \
                $BAM;
        echo "Indexing sorted BAM $SORTED_BAM";
        sambamba index \
                --nthreads $NUMBER_PROCESSORS \
                --show-progress \
                $BAM;

}

download_and_index_reference_genome;
align_fastq_pairs $NORMAL_FASTQ_DIR $NORMAL_SAMPLE_NAME;
align_fastq_pairs $TUMOR_FASTQ_DIR TUMOR_SAMPLE_NAME;
process_alignments "$NORMAL_SAMPLE_NAME.bam";
process_alignments "$TUMOR_SAMPLE_NAME.bam";

