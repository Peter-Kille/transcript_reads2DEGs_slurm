#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=100

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo "\$SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "\$SLURM_NTASKS=${SLURM_NTASKS}"
echo "\$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}"
echo "\$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}"
echo "\$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}"
echo "\$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}"

# read in sample and read names and move and rename samples into workdir
tail -n +2 ${pipedir}/"${metadata}" | while IFS=$'\t' read -r sample read1 read2 rest; do

echo ${sample}
echo ${read1}
echo ${read2}
cp ${sourcedir}/${read1} ${rawdir}/${sample}_1.fastq.gz
cp ${sourcedir}/${read2} ${rawdir}/${sample}_2.fastq.gz

done

awk -F'\t' 'BEGIN { OFS="\t" }
NR==1 { print $1, "Files", substr($0, index($0,$2)) }
NR>1 { print $1, $1 ".markdup.featurecount", substr($0, index($0,$2)) }
' ${pipedir}/${metadata} > ${pipedir}/${metadata}.fc

