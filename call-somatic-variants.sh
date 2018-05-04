#!/bin/bash
set -e

# machine configuration
NUMBER_PROCESSORS=32
MEMORY_LIMIT=60GB

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

    if [ $# -ne 1 ] ; then
        echo "Function 'run' expected 1 argument but got $#";
        exit 1;
    fi

    local COMMAND=$1;

    echo $COMMAND;
    eval "time $COMMAND"
}

function run_unless_exists() {
    # Takes 3 arguments:
    #   1) Description of what work might get executed
    #   2) File name, which if it exists and is not empty, should prevent
    #       anything further from being executed.
    #   3) Command to run.

    if [ $# -ne 3 ] ; then
        echo "Function 'run_unless_exists' expected 3 arguments but got $#";
        exit 1;
    fi

    # run a command unless its output file has already been created
    local DESCRIPTION=$1;
    local OUTPUT_FILE=$2;
    local COMMAND=$3;

    if [ -s $OUTPUT_FILE ]; then
        echo ">> Skipping '$DESCRIPTION' because $OUTPUT_FILE already exists";
    else
        echo ">> Running '$DESCRIPTION' to create $OUTPUT_FILE";
        run "$COMMAND";
    fi
}

function download_and_index_reference_genome() {
    if [ $# -ne 0 ] ; then
        echo "Function 'download_and_index_reference_genome' expected 0 arguments but got $#";
        exit 1;
    fi

    run_unless_exists "Download reference" $REFERENCE_FASTA_PATH \
        "wget $REFERENCE_FASTA_SOURCE && gunzip hs37d5.fa.gz && mv hs37d5.fa $REFERENCE_FASTA_PATH";
    run_unless_exists "Creating index for $REFERENCE_FASTA_PATH" $REFERENCE_INDEX_PATH \
        "bwa index $REFERENCE_FASTA_PATH";
}


function align_fastq_pairs() {
    # Align every FASTQ pair into multiple BAM files
    local FASTQ_DIR=$1;
    local FASTQ_PREFIX=$2;

    if [ $# -ne 2 ] ; then
        echo "Function 'align_fastq_pairs' expected 2 arguments but got $#";
        exit 1;
    fi

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
            local BAM="$READ_GROUP.bam"
            local READ_GROUP_TAG="'@RG\tID:$READ_GROUP\tSM:$FASTQ_PREFIX\tLB:$FASTQ_PREFIX\tPL:ILLUMINA'"
            run_unless_exists "Align $R1_fastq and $R2_fastq" $BAM \
                "bwa mem -M \
                        -t $NUMBER_PROCESSORS \
                        -R $READ_GROUP_TAG \
                        $REFERENCE_FASTA_PATH \
                        $R1_fastq \
                        $R2_fastq \
                        | samtools view -S -b -@$NUMBER_PROCESSORS -o $BAM -";
    done
}

function process_alignments() {
    # Runs the following pipeline steps:
    #       - sort BAM
    #       - index BAM
    #       - mark duplicates
    # Input: sample bam (expected to exist $SAMPLE_NAME.bam)
    # Output: generates SAMPLE_NAME.final.bam

    if [ $# -ne 1 ] ; then
        echo "Function 'process_alignments' expected 1 argument but got $#";
        exit 1;
    fi

    local UNSORTED_BAM_PREFIX=$1;

    # sort and index every BAM
    # for any name X.bam, create X.sorted.bam
    # exclude file names which contain sequence '.sorted.'
    for UNSORTED_BAM in $UNSORTED_BAM_PREFIX*.bam ; do
        case $UNSORTED_BAM in
            *.sorted.*)
                echo "Skipping $UNSORTED_BAM since it contains '.sorted.'";
                continue;;
            *.merged.*)
                echo "Skipping $UNSORTED_BAM since it contains '.merged.'";
                continue;;
            *.final.*)
                echo "Skipping $UNSORTED_BAM since it contains '.final.'";
                continue;;
            *)
                local SORTED_BAM=`echo $UNSORTED_BAM | sed -e 's/\.bam/\.sorted\.bam/g'`
                run_unless_exists "Sorting $UNSORTED_BAM" $SORTED_BAM \
                    "sambamba sort \
                            --memory-limit $MEMORY_LIMIT \
                            --show-progress \
                            --nthreads $NUMBER_PROCESSORS \
                            --out $SORTED_BAM \
                            $UNSORTED_BAM";

                local SORTED_BAM_INDEX="$SORTED_BAM.bai"
                run_unless_exists "Indexing sorted BAM $SORTED_BAM" $SORTED_BAM_INDEX \
                    "sambamba index \
                        --nthreads $NUMBER_PROCESSORS \
                        --show-progress \
                        $SORTED_BAM";
        esac
    done
    local MERGED_BAM="$UNSORTED_BAM_PREFIX.merged.bam"
    run_unless_exists "Merging lanes" $MERGED_BAM \
        "sambamba merge \
            --nthreads $NUMBER_PROCESSORS \
            --show-progress \
            $MERGED_BAM \
            $UNSORTED_BAM_PREFIX*.sorted.bam";
    local FINAL_BAM="$UNSORTED_BAM_PREFIX.final.bam";

    # for larger WGS, need to have both larger overflow
    # list and hash table to avoid hitting too many open
    # files
    run_unless_exists "Marking duplicates" $FINAL_BAM \
        "sambamba markdup \
                --nthreads $NUMBER_PROCESSORS \
                --show-progress \
                --overflow-list-size 1000000 \
                --hash-table-size 4194304 \
                $MERGED_BAM \
                $FINAL_BAM";
    local FINAL_BAM_INDEX="$FINAL_BAM.bai"
    run_unless_exists "Indexing final BAM" $FINAL_BAM_INDEX \
        "sambamba index \
                --nthreads $NUMBER_PROCESSORS \
                --show-progress \
                $FINAL_BAM";
}


function call_somatic_variants() {
    if [ $# -ne 2 ] ; then
        echo "Function 'call_somatic_variants' expected 2 argument but got $#";
        exit 1;
    fi

    local NORMAL_FASTQ_PREFIX=$1;
    local TUMOR_FASTQ_PREFIX=$2;
    local NORMAL_BAM="$NORMAL_FASTQ_PREFIX.final.bam";
    local TUMOR_BAM="$TUMOR_FASTQ_PREFIX.final.bam";

    # Strelka2 expects a .fai file associated with the reference
    run_unless_exists "Indexing reference FASTA" "$REFERENCE_FASTA_PATH.fai" \
        "samtools faidx $REFERENCE_FASTA_PATH";

    local VCF_PREFIX="$NORMAL_FASTQ_PREFIX.$TUMOR_FASTQ_PREFIX"
    local STRELKA_DIR="Strelka.$VCF_PREFIX"

    run_unless_exists "Generating Strelka2 configuration" "$STRELKA_DIR/runWorkflow.py" \
        "mkdir -p $STRELKA_DIR \
         && configureStrelkaSomaticWorkflow.py \
            --normalBam $NORMAL_BAM \
            --tumorBam $TUMOR_BAM \
            --referenceFasta $REFERENCE_FASTA_PATH \
            --runDir $STRELKA_DIR";

    # execution on a single local machine with multiple threads
    run_unless_exists "Calling somatic variants" "$STRELKA_DIR/results/stats/runStats.tsv" \
        "cd $STRELKA_DIR \
         && python runWorkflow.py \
                -m local \
                -j $NUMBER_PROCESSORS \
                -g $MEMORY_LIMIT \
         && cd .."

    local SNV_VCF="$VCF_PREFIX.snvs.vcf"
    run_unless_exists "Decompressing and renaming SNV VCF" $SNV_VCF \
        "zcat $STRELKA_DIR/results/variants/somatic.snvs.vcf.gz > $SNV_VCF";
    local INDEL_VCF="$VCF_PREFIX.indels.vcf"
    run_unless_exists "Decompressing and renaming indel VCF" $INDEL_VCF \
        "zcat $STRELKA_DIR/results/variants/somatic.indels.vcf.gz > $INDEL_VCF";
    echo "Summary:";
    echo "  Number of passing SNVs: `cat $SNV_VCF | grep PASS | wc -l`"
    echo "  Number of passing indels: `cat $INDEL_VCF | grep PASS | wc -l`"
}

download_and_index_reference_genome;
align_fastq_pairs $NORMAL_FASTQ_DIR $NORMAL_FASTQ_PREFIX;
align_fastq_pairs $TUMOR_FASTQ_DIR $TUMOR_FASTQ_PREFIX;
process_alignments $NORMAL_FASTQ_PREFIX;
process_alignments $TUMOR_FASTQ_PREFIX;
call_somatic_variants $NORMAL_FASTQ_PREFIX $TUMOR_FASTQ_PREFIX;
