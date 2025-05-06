#!/bin/bash

# Source config script
source arguments_param
source config.parameters

# Step 1: Rename and count file
# -- Copy raw data files from sourcedir to rawdir.
# CORE PARAMETERS: sourcedir, rawdir
# INPUT: 
# WORK: rawdir
# OUTPUT: null
# PROCESS - File transfer
sbatch -d singleton --error="${log}/1-preprocess_%J.err" --output="${log}/1-preprocess_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/1-preprocess.sh"

samples=$( tail -n +2 ${pipedir}/${metadata} | cut -f1 )

export samples
export sample_array=($samples)
sample_number=${#sample_array[@]}
sample_number=$(($sample_number - 1))

echo $samples
echo $sample_number

# Check for zero files
if [[ "$sample_number" -eq 0 ]]; then
    echo "Error: No files found in ${metadata}"
    exit 1
fi

exit

# Step 2: QC - fastqc, fastp, fastqc
# -- Run FastQC on raw data to assess data quality before trimming.
# CORE PARAMETERS: modules, rawdir, qcdir, log
# INPUT: rawdir
# WORK: qcdir
# OUTPUT: null
sbatch -d singleton --error="${log}/2A-rawqc_%J.err" --output="${log}/2A-rawqc_%J.out" --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2A-fastqc_array.sh"

sbatch -d singleton --error="${log}/2B-fastp_%J.err" --output="${log}/2B-fastp_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2B-fastp_array.sh"  

sbatch -d singleton --error="${log}/2C-trimqc_%J.err" --output="${log}/2C-trimqc_%J.out" --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2C-fastqc-trim.sh"

# Step 3: star genome indexing and mapping 
# DESCRIPTION: genome indexing and mapping
# CORE PARAMETERS: modules, rawdir, trimdir, stardir, genomedir
# INPUT: samples, trimdir, genomedir
# WORK: genomedir, stardir
# OUTPUT: null
sbatch -d singleton --error="${log}/3A_star_index_%J.err" --output="${log}/3A_star_index_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/3A-star_index_genome.sh"

sbatch -d singleton --error="${log}/3B-star_map_%J.err" --output="${log}/3B-star_map_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/3B-star_mapping_array.sh"

# Step 4: Mark duplicates
# DESCRIPTION: Mark Duplicates
# CORE PARAMETERS: modules, trimdir, log, markdir, samples
# INPUT: trimdir, samples
# WORK: trimdir, markdir
# OUTPUT: null

sbatch -d singleton --error="${log}/4-markdup_%J.err" --output="${log}/4-markdup_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/4-markdup_array.sh"

# Step 5: featurecount
# DESCRIPTION: feature count with subread
# CORE PARAMETERS: modules, trimdir, log, samples, fcdir, markdir
# INPUT: samples, trimdir, markdir
# WORK: markdir, fcdir 
# OUTPUT: featurecounts

sbatch -d singleton --error="${log}/5-featurecount_%J.err" --output="${log}/5-featurecount_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/5-featurecount_array.sh"

# Step X: MultiQC report
# -- Generate a MultiQC report to summarize the results of all previous steps.
# CORE PARAMETERS: modules, workdir, multiqc
# INPUT: workdir
# WORK: multiqc
# OUTPUT: multiqc
# PROCESS - multiqc
sbatch -d singleton --error="${log}/multiqc_%J.err" --output="${log}/multiqc_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/X-multiqc.sh"

# Step Z: template
# DESCRIPTION:.
# CORE PARAMETERS: modules, scripts rawdir, trimdir, log
# INPUT: rawdir, trimdir
# WORK: metaphlandir, metaphladbdir
# OUTPUT: null
# PROCESS - map raw reads against database for biodiversoty classification
# export
#sbatch -d singleton --error="${log}/8A-metaphlandb_%J.err" --output="${log}/8A-metaphlandb_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/8A-metaphlandb_download.sh"

#sbatch -d singleton --error="${log}/8B-metaphlan_%J.err" --output="${log}/8B-metaphlan_%J.out" --array="1-${sample_number}%8" --job-name=${NAME} --partition=${PART} "${moduledir}/8B-metaphlan_array.sh"

#sbatch -d singleton --error="${log}/8C-metaphlan_post_%J.err" --output="${log}/8C-metaphlan_post_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/8C-metaphan_post.sh"
