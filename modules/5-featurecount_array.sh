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

# Load some modules
module load ${subread_module}

sample_array=($samples)
base=${sample_array[$SLURM_ARRAY_TASK_ID]}

featureCounts \
        -T 4 -p -F GTF -t exon -g gene_id \
	-a ${genomedir}/${annot} \
	-o ${fcdir}/${base}.markdup.featurecount \
	${markdir}/${base}.markdup.bam

featureCounts \
	-T 4 -p -F GTF -t exon -g gene_id \
	-a ${genomedir}/${annot} \
	-o ${fcdir}/${base}.rmdup.featurecount \
	${markdir}/${base}.rmdup.bam

cp -R ${fcdir} ${outdir}/

