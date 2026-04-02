#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o assembly.test.log
#SBATCH --account=nsyamat8807
#SBATCH --partition=silver

module load biological/perl_5.40
module load biological/samtools_1.23
module load biological/java

export PROJ_DIR=/export/home/bioclass/nsyamat8807/labnine
cd $PROJ_DIR
export SRR=SRR5324768


# worked
java -jar /export/share/software/biological/picard/picard.jar \
	CreateSequenceDictionary \
    REFERENCE=genome/Thermus_thermophilus_TTHNAR1.fa \
    OUTPUT=genome/Thermus_thermophilus_TTHNAR1.dict


# worked
 /export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2-build \
	genome/Thermus_thermophilus_TTHNAR1.fa \
	genome/Thermus_thermophilus_TTHNAR1


# run bowtie2 and samtools separately for testing. bowtie2 failed.
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 -x \
		genome/Thermus_thermophilus_TTHNAR1 \
        -1 fastq/${SRR}_pass_1.fastq.gz \
        -2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
        --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
        > alignment/${SRR}.sam 


# worked
 samtools view -hb alignment/${SRR}.sam | samtools sort -l 5 -o alignment/${SRR}.bam




samtools consensus -f fasta alignment/${SRR}.bam -o _consensus.fasta