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

echo "Simple somatic variant calling pipeline";
echo "=======";

if [ $# -ne 4 ] ; then
    echo "Wrong number of arguments ($#)";
    echo "Expected arguments:";
    echo "      (1) directory containing normal FASTQ files";
    echo "      (2) common prefix in names of all normal FASTQ files";
    echo "      (3) directory containing tumor FASTQ files";
    echo "      (4) common prefix in names of all tumor FASTQ files";
    echo "----"
    echo "Example:";
    echo "      ./call-somatic-variants.sh . normal . tumor";
    exit 1;
else
    NORMAL_FASTQ_DIR=$1;
    NORMAL_FASTQ_PREFIX=$2;
    TUMOR_FASTQ_DIR=$3;
    TUMOR_FASTQ_PREFIX=$4;
fi

echo "Quick & Dirty Somatic Variant Calling Pipeline";
echo "---";
echo "Normal FASTQ location: $NORMAL_FASTQ_DIR/$NORMAL_FASTQ_PREFIX*.fastq";
echo "Tumor FASTQ location: $TUMOR_FASTQ_DIR/$TUMOR_FASTQ_PREFIX*.fastq";
echo "---";


function run() {
        # print a command before running it wrapped in a 'time' command
        local COMMAND=$1;
        echo $COMMAND;
        eval "time $COMMAND"
}

function download_and_index_reference_genome() {
        echo "-- download_and_index_reference_genome";
        if [ ! -e $REFERENCE_FASTA_PATH ]; then
                echo "Couldn't find reference file $REFERENCE_FASTA_PATH";
                echo "Downloading from $REFERENCE_FASTA_SOURCE..."
                run "wget $REFERENCE_FASTA_SOURCE";
                echo "Decompressing downloaded reference genome..."
                run "gunzip hs37d5.fa.gz";
                run "mv hs37d5.fa $REFERENCE_FASTA_PATH";
        else
                echo "Using reference: $REFERENCE_FASTA_PATH"
        fi;
        if [ ! -s $REFERENCE_INDEX_PATH ]; then
                echo "Creating index for $REFERENCE_FASTA_PATH"
                run "bwa index $REFERENCE_FASTA_PATH";
        fi;
}


function align_fastq_pairs() {
        # 1) Align every FASTQ pair into multiple BAM files
        # 2) Merge them into single BAM file per sample
        local FASTQ_DIR=$1;
        local FASTQ_PREFIX=$2;
        echo "-- align_fastq_pairs";
        echo "  FASTQ_DIR=$FASTQ_DIR";
        echo "  FASTQ_PREFIX=$FASTQ_PREFIX";

        # check to make sure that all arguments are non-empty
        if [[ -z $FASTQ_DIR ]] ; then
                echo "Missing first argument (FASTQ_DIR)";
                exit 1;
        fi
        if [[ -z $FASTQ_PREFIX ]] ; then
                echo "Missing second argument (FASTQ_PREFIX)";
                exit 1;
        fi

        for R1_fastq in $FASTQ_DIR/$FASTQ_PREFIX*.R1.fastq.gz ; do
                R2_fastq=`echo $R1_fastq | sed -e 's/\.R1\./\.R2\./g'`
                if [ ! -e $R2_fastq ]; then
                        echo "Couldn't find R2 ($R2_fastq) corresponding to $R1_fastq"
                        exit 1;
                fi;
                echo "R1: $R1_fastq";
                echo "R2: $R2_fastq";
                # make a local file name for the BAM we're going to generate from each FASTQ pair
                local READ_GROUP=`basename $R1_fastq | sed -e 's/\.R1\.fastq\.gz//g'`
                local SAMPLE=
                local SAM="$READ_GROUP.sam"
                local BAM="$READ_GROUP.bam"
                # test if $BAM exists and is non-empty
                if [ ! -s $BAM ]; then
                        echo "Generating BAM file $BAM";
                        run "bwa mem -M \
                                -t $NUMBER_PROCESSORS \
                                -R '@RG\tID:$READ_GROUP\tSM:$FASTQ_PREFIX\tLB:$FASTQ_PREFIX\tPL:ILLUMINA' \
                                $REFERENCE_FASTA_PATH \
                                $R1_fastq $R2_fastq > $SAM";
                        run "samtools view -S -b $SAM \
                                -@ $NUMBER_PROCESSORS \
                                -o $BAM";
                        # once we've successfully create BAM file, clean up SAM
                        run "rm $SAM";
                        echo "---";
                else
                        echo "BAM file $BAM already exists";
                fi;
        done;
        echo "Merging all BAM files for $FASTQ_PREFIX into single alignment file";
        bamtools merge -in $FASTQ_DIR/$FASTQ_PREFIX*.bam  -out $FASTQ_PREFIX.bam;
}

function process_alignments() {
        # Runs the following pipeline steps:
        #       - sort BAM
        #       - index BAM
        #       - mark duplicates
        # Input: sample bam (expected to exist $SAMPLE_NAME.bam)
        # Output: generates SAMPLE_NAME.final.bam

        # sort and index BAM file
        local UNSORTED_BAM=$1;
        local SORTED_BAM=`echo $UNSORTED_BAM | sed -e 's/\.bam/\.sorted\.bam/g'`
        echo "-- process_alignments";
        echo "  UNSORTED_BAM: $UNSORTED_BAM";
        echo "  SORTED_BAM: $SORTED_BAM";

        echo "Sorting $BAM to generate $SORTED_BAM";
        run "sambamba sort \
                --memory-limit $MEMORY_LIMIT \
                --show-progress \
                --nthreads $NUMBER_PROCESSORS \
                --out $SORTED_BAM \
                $BAM";
        echo "Indexing sorted BAM $SORTED_BAM";
        run "sambamba index \
                --nthreads $NUMBER_PROCESSORS \
                --show-progress \
                $BAM";
        # strip off everything after the first '.' in the file name
        local BASE_FILENAME=`basename $BAM | cut -f 1 -d '.'`;
        local FINAL_BAM="$BASE_FILENAME.final.bam";
        echo "Marking duplicates to generate $FINAL_BAM";
        run "sambamba markdup \
                --nthreads $NUMBER_PROCESSORS \
                --show-progress \
                $SORTED_BAM \
                $FINAL_BAM";
}


function call_somatic_variants() {
        local NORMAL_BAM=$1;
        local TUMOR_BAM=$2;
        echo "-- call_somatic_variants";
        echo "  NORMAL_BAM: $NORMAL_BAM";
        echo "  TUMOR_BAM: $TUMOR_BAM";
        echo "Generating Strelka2 configuration";
        run "configureStrelkaSomaticWorkflow.py \
                --normalBam $NORMAL_BAM \
                --tumorBam $TUMOR_BAM \
                --referenceFasta $REFERENCE_FASTA_PATH \
                --runDir .";
        echo "Running Strelka2";
        # execution on a single local machine with 20 parallel jobs
        run "./runWorkflow.py -m local -j $NUMBER_PROCESSORS";
}

download_and_index_reference_genome;
align_fastq_pairs $NORMAL_FASTQ_DIR $NORMAL_FASTQ_PREFIX;
align_fastq_pairs $TUMOR_FASTQ_DIR $TUMOR_FASTQ_PREFIX;
process_alignments "$NORMAL_FASTQ_PREFIX.bam";
process_alignments "$TUMOR_FASTQ_PREFIX.bam";
call_somatic_variants "$NORMAL_FASTQ_PREFIX.final.bam" "$TUMOR_FASTQ_PREFIX.final.bam";
