#!/bin/bash
#SBATCH -n 1                # Number of tasks to run each step
#SBATCH -N 1                # Number of nodes requested
#SBATCH --cpus-per-task=12     # Number of cpu cores per task
#SBATCH -t 24:00:00            # Runtime in minutes (36 hours max)
#SBATCH --qos=medium         # The QoS to submit the job.
#SBATCH --mem=2G           # Memory per cpu in G (see also --mem-per-cpu)
#SBATCH -o mapping.out  # Standard output goes to this file
#SBATCH -e mapping.err  # Standard error goes to this file

# Import modules
module load bowtie2
module load samtools
module load fastqc

## PATHS and variables for mapping
GenDB="./Ref_genome/ecoli"
FQdir="./seq_files/"
QCoutdir="/home/romonmar/alvaro"
thread=12
NAMES="./filenames.txt" # A tuned file with the names of the files stopping in the last _
R1="_1.fq"
R2="_2.fq"
OUTPUT_DIR="./bam_files/"
BAMEND="_sort.bam"

# List all fq.gz files and for each one start unzipping them
for file in $FQdir*
do
	# The gzip can't be parallelized, only uses 1 cpu and 1 node
	srun -n 1 -c 1 gzip -d $file
	wait
	# Take de $file variable, eliminate .gz extension and use fastq
	# for an html report using 12 cpus
	file=${file:0:-3}
	srun -n 1 -c 12 fastqc -t $thread $file --outdir $QCoutdir
	wait
done


# List all samples names and start mapping with bowtie2
for file in $(cat $NAMES)
do
	# bowtie2 mapping
	srun -n 1 -c 12 bowtie2 -x $GenDB -1 $FQdir$file$R1 -2 $FQdir$file$R2 --no-unal -p $thread -S $OUTPUT_DIR$file.sam
	wait
	# From sam to bam
	srun -n 1 -c 1 samtools view -bSq 30 $OUTPUT_DIR$file.sam > $OUTPUT_DIR$file.bam
	wait
	rm $OUTPUT_DIR$file.sam # Free disk space
	# Construct sorted bam
	srun -n 1 -c 1 samtools sort $OUTPUT_DIR$file.bam -o $OUTPUT_DIR$file$BAMEND
        wait
	rm $OUTPUT_DIR$file.bam
done
