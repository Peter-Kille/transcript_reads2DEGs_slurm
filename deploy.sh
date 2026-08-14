#!/bin/bash

###### Source config script ######
source config/arguments_param
source config/folders
source config/programs

###### Step 1: Rename and count file
sbatch -d singleton --error="${log}/1-preprocess_%J.err" --output="${log}/1-preprocess_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/1-preprocess.sh"

samples=$( tail -n +2 ${sourcedir}/${metadata} | cut -f1,1 )

export samples
echo $samples
export sample_array=($samples)
sample_number=${#sample_array[@]}
echo $sample_number
sample_number=$(($sample_number - 1))

####### Check for zero files
if [[ "$sample_number" -eq 0 ]]; then
    echo "Error: No files found in ${metadata}"
    exit 1
fi

###### Step 2: QC - fastqc, fastp, fastqc
sbatch -d singleton --error="${log}/2A-rawqc_%J.err" --output="${log}/2A-rawqc_%J.out" --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2A-fastqc_array.sh"

sbatch -d singleton --error="${log}/2B-fastp_%J.err" --output="${log}/2B-fastp_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2B-fastp_array.sh"  

sbatch -d singleton --error="${log}/2C-trimqc_%J.err" --output="${log}/2C-trimqc_%J.out" --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/2C-fastqc-trim.sh"

###### Step 3: star genome indexing and mapping 

sbatch -d singleton --error="${log}/3A-star_index_%J.err" --output="${log}/3A-star_index_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/3A-star_index_genome.sh"

sbatch -d singleton --error="${log}/3B-star_map_%J.err" --output="${log}/3B-star_map_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/3B-star_mapping_array.sh"

###### Step 4: MapQC

sbatch -d singleton --error="${log}/4-mapqc_%J.err" --output="${log}/4-mapqc_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/4-mapqc_array.sh"

###### Step 5: featurecount
sbatch -d singleton --error="${log}/5-featurecount_%J.err" --output="${log}/5-featurecount_%J.out" --"array=0-${sample_number}%20" --job-name=${NAME} --partition=${PART} "${moduledir}/5-featurecount_array.sh"

###### Step 6: SARTools

sbatch -d singleton --error="${log}/6-sartools_%J.err" --output="${log}/6-sartools_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/6-SARTools.sh"

###### Step X: MultiQC report
sbatch -d singleton --error="${log}/multiqc_%J.err" --output="${log}/multiqc_%J.out" --job-name=${NAME} --partition=${PART} "${moduledir}/X-multiqc.sh"
