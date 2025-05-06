#!/bin/bash
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=8      #
#SBATCH --mem-per-cpu=2000     # in megabytes, unless unit explicitly stated

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

#load some modules
module load ${picard_module}
module load ${samtools_module}

sample_array=($samples)
base=${sample_array[$SLURM_ARRAY_TASK_ID]}

samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $stardir/${base}.sorted.bam $stardir/${base}-unsort.Aligned.out.bam
samtools index $stardir/${base}.sorted.bam

##  MARK DUPLICATES  ##
picard MarkDuplicates I=$stardir/${base}.sorted.bam O=$markdir/${base}.markdup.bam M=$markdir/${base}.metrics.markdup.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

## REMOVE DUPLICATES ##
picard MarkDuplicates I=${stardir}/${base}.sorted.bam O=${markdir}/${base}.rmdup.bam M=${markdir}/${base}.metrics.rmdup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

