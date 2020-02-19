#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH -o /data/biohub/20190717_Lardelli_RNASeq_Larvae/slurm/%x_%j.out
#SBATCH -e /data/biohub/20190717_Lardelli_RNASeq_Larvae/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stephen.pederson@adelaide.edu.au

## Cores
CORES=16

## Modules
module load FastQC/0.11.7
module load STAR/2.5.3a-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load AdapterRemoval/2.2.1-foss-2016b
module load GCC/5.4.0-2.26
module load Subread

## Genomic Data Files
REFS=/data/biorefs/reference_genomes/ensembl-release-96/danio_rerio
if [ ! -d "${REFS}" ]; then
  echo -e "Couldn't find ${REFS}"
  exit 1
fi
echo -e "Found ${REFS}"

GTF=${REFS}/Danio_rerio.GRCz11.96.chr.gtf.gz
if [ ! -f "${GTF}" ]; then
  echo -e "Couldn't find ${GTF}"
  exit 1
fi
echo -e "Found ${GTF}"


## Directories
PROJROOT=/data/biohub/20190717_Lardelli_RNASeq_Larvae
if [ ! -d "${PROJROOT}" ]; then
  echo -e "Couldn't find ${PROJROOT}"
  exit 1
fi
echo -e "Found ${PROJROOT}"

## Directories for Initial FastQC
RAWDATA=${PROJROOT}/0_rawData
mkdir -p ${RAWDATA}/FastQC
if [ ! -d "${RAWDATA}/fastq" ]; then
  echo -e "Couldn't find ${RAWDATA}/fastq"
  exit 1
fi
echo -e "Found ${RAWDATA}/fastq"

## Setup for Trimmed data
TRIMDATA=${PROJROOT}/1_trimmedData
mkdir -p ${TRIMDATA}/fastq
mkdir -p ${TRIMDATA}/FastQC
mkdir -p ${TRIMDATA}/log

## Setup for genome alignment
ALIGNDATA=${PROJROOT}/2_alignedData
mkdir -p ${ALIGNDATA}/log
mkdir -p ${ALIGNDATA}/bam
mkdir -p ${ALIGNDATA}/FastQC
mkdir -p ${ALIGNDATA}/featureCounts


##--------------------------------------------------------------------------------------------##
## FastQC on the raw data
##--------------------------------------------------------------------------------------------##

fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/*.fq.gz

##--------------------------------------------------------------------------------------------##
## Trimming the Merged data
##--------------------------------------------------------------------------------------------##

for R1 in ${RAWDATA}/fastq/*1.fq.gz
  do

    echo -e "Currently working on ${R1}"

    # Now create the output filenames
    out1=${TRIMDATA}/fastq/$(basename $R1)
    BNAME=${TRIMDATA}/fastq/$(basename ${R1%_R1.fq.gz})
    echo -e "Output file 1 will be ${out1}"
    echo -e "Trimming:\t${BNAME}"

    #Trim
    AdapterRemoval \
      --gzip \
      --trimns \
      --trimqualities \
      --minquality 30 \
      --minlength 35 \
      --threads ${CORES} \
      --basename ${BNAME} \
      --output1 ${out1} \
      --file1 ${R1}

  done

## Move the log files into their own folder
mv ${TRIMDATA}/fastq/*settings ${TRIMDATA}/log

## Run FastQC
fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastq/*.fq.gz


##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to the genome
##--------------------------------------------------------------------------------------------##

## Aligning, filtering and sorting
for R1 in ${TRIMDATA}/fastq/*1.fq.gz
  do

  BNAME=$(basename ${R1%_R1.fq.gz})
  echo -e "STAR will align:\t${R1}"

    STAR \
        --runThreadN ${CORES} \
        --genomeDir ${REFS}/star \
        --readFilesIn ${R1} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix ${ALIGNDATA}/bam/${BNAME} \
        --outSAMtype BAM SortedByCoordinate

  done

## Move the log files into their own folder
mv ${ALIGNDATA}/bam/*out ${ALIGNDATA}/log
mv ${ALIGNDATA}/bam/*tab ${ALIGNDATA}/log

## Fastqc and indexing
for BAM in ${ALIGNDATA}/bam/*.bam
 do
   fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATA}/FastQC --noextract ${BAM}
   samtools index ${BAM}
 done

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
sampleList=`find ${ALIGNDATA}/bam -name "*out.bam" | tr '\n' ' '`

## featureCounts needs an extracted gtf
zcat ${GTF} > temp.gtf

## Running featureCounts on the sorted bam files
featureCounts -Q 10 \
  -s 1 \
  -T ${CORES} \
  --fracOverlap 1 \
  -a temp.gtf \
  -o ${ALIGNDATA}/featureCounts/counts.out ${sampleList}

## Remove the temp gtf
rm temp.gtf

## Storing the output in a single file
cut -f1,7- ${ALIGNDATA}/featureCounts/counts.out | sed 1d > ${ALIGNDATA}/featureCounts/genes.out
