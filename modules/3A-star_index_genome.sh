#!/bin/bash
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=32      #
#SBATCH --mem=128000     # in megabytes, unless unit explicitly stated

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
module load ${star_module}

cp ${pipedir}/genome/* ${genomedir}/

## Change --sjdbOverhang to length of your sequence data (one read of your paired data) minus 1

STAR 	--runThreadN ${SLURM_CPUS_PER_TASK} \
        --limitGenomeGenerateRAM 128000000000 \
	--runMode genomeGenerate \
	--genomeDir  ${genomedir} \
	--genomeFastaFiles ${genomedir}/${genome} \
	--sjdbGTFfile ${genomedir}/${annot} \
	--sjdbOverhang ${len}
