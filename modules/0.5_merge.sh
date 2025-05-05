#!/bin/bash
#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=16      #   
#SBATCH --mem-per-cpu=500     # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err         # redirect stderr to this file
#SBATCH --output=%J.out        # redirect stdout to this file
##SBATCH --mail-user=[insert email address]@Cardiff.ac.uk  # email address used for event notification
##SBATCH --mail-type=end                                   # email on job end
##SBATCH --mail-type=fail                                  # email on job failure

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

# Write jobscript to output file (good for reproducibility)
cat $0

workdir=$(pwd)
rawdata=${workdir}/rawdata
#create directory for qc output
mkdir ${workdir}/merge
mergedir=${workdir}/merge

for f in ${rawdata}/*_L001_R1_001.fastq.gz
do

R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_L001_R1_001//')

#pigz -p ${SLURM_CPUS_PER_TASK} -9 ${rawdata}/${base}_L001_R1_001.fastq
#pigz -p ${SLURM_CPUS_PER_TASK} -9 ${rawdata}/${base}_L002_R1_001.fastq
#pigz -p ${SLURM_CPUS_PER_TASK} -9 ${rawdata}/${base}_L001_R2_001.fastq
#pigz -p ${SLURM_CPUS_PER_TASK} -9 ${rawdata}/${base}_L002_R2_001.fastq

cat ${rawdata}/${base}_L001_R1_001.fastq.gz ${rawdata}/${base}_L002_R1_001.fastq.gz > ${mergedir}/${base}_R1.fastq.gz
cat ${rawdata}/${base}_L001_R2_001.fastq.gz ${rawdata}/${base}_L002_R2_001.fastq.gz > ${mergedir}/${base}_R2.fastq.gz

done
