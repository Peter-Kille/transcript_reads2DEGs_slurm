#!/bin/bash

# deploy.sh plays two roles (see PIPELINE_GUIDE.md in metagenomic_slurm):
#   launcher - validates options on the login node, submits itself as a small
#              Slurm "driver" job and exits
#   driver   - runs on a compute node (PIPELINE_DRIVER=1), submits one stage at a
#              time with `sbatch --wait` and records a done marker per stage

# The driver starts with --mem=1G; stop that leaking into child sbatch calls
if [[ -n "${PIPELINE_DRIVER}" ]]; then
    unset SLURM_MEM_PER_NODE SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU
fi

###### Source config script ######
source config/arguments_param
source config/folders
source config/programs

###### Record the command line (once per launch, not again in the driver) ######
if [[ -z "${PIPELINE_DRIVER}" ]]; then {
    echo "# $(date '+%Y-%m-%d %H:%M:%S') user=${USER} host=$(hostname)"
    echo "# cwd: $(pwd -P)"
    echo "${CMDLINE}"
    echo ""
} >> "${log}/0_execution_commandline.txt"; fi

###### Launcher: hand over to the driver job and exit ######
if [[ -z "${PIPELINE_DRIVER}" && "${foreground}" != "yes" ]]; then
    driver_id=$(sbatch --parsable \
        --job-name="${NAME}_driver" --partition=${driver_part} \
        --time="${driver_time}" --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=1G \
        --output="${log}/0_driver_%J.out" --error="${log}/0_driver_%J.err" \
        --export=ALL,PIPELINE_DRIVER=1 \
        --wrap "cd '$(pwd -P)' && bash ${CMDLINE}")
    echo "Driver job submitted: ${driver_id}"
    echo "Follow progress: tail -f ${log}/0_driver_${driver_id}.out"
    exit 0
fi

###### run(): submit one stage, wait for it, write a done marker ######
# Usage: run STEP_NAME <sbatch options...> <module script>
run() {
    step_name="$1"
    shift

    marker="${workdir}/.done_${step_name}"

    # Already finished in an earlier run? Skip it.
    if [ -e "$marker" ]; then
        echo "$(date '+%F %T') SKIP   ${step_name} (done marker found)"
        return 0
    fi

    echo "$(date '+%F %T') SUBMIT ${step_name}"
    sbatch --wait "$@"
    exit_code=$?

    if [ "$exit_code" -eq 0 ]; then
        touch "$marker"
        echo "$(date '+%F %T') DONE   ${step_name}"
    else
        echo "$(date '+%F %T') FAILED ${step_name} - stopping the pipeline" >&2
        exit 1
    fi
}

###### Step 1: Rename and count file
run 1-preprocess \
    --error="${log}/1-preprocess_%J.err" --output="${log}/1-preprocess_%J.out" \
    --job-name=${NAME} --partition=${PART} "${moduledir}/1-preprocess.sh"

# Count samples only after step 1 has finished (run waits, so this is safe)
samples=$( tail -n +2 ${sourcedir}/${metadata} | cut -f1,1 )

export samples
echo $samples
export sample_array=($samples)
sample_number=${#sample_array[@]}
echo $sample_number

####### Check for zero files
if [[ "$sample_number" -eq 0 ]]; then
    echo "Error: No files found in ${metadata}" >&2
    exit 1
fi
sample_number=$(($sample_number - 1))

###### Steps 2 and 3A: QC / trimming and STAR genome index run in parallel
prep_pids=()

(
    run 2A-fastqc_array \
        --error="${log}/2A-rawqc_%J.err" --output="${log}/2A-rawqc_%J.out" \
        --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} \
        "${moduledir}/2A-fastqc_array.sh"

    run 2B-fastp_array \
        --error="${log}/2B-fastp_%J.err" --output="${log}/2B-fastp_%J.out" \
        --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} \
        "${moduledir}/2B-fastp_array.sh"

    run 2C-fastqc-trim \
        --error="${log}/2C-trimqc_%J.err" --output="${log}/2C-trimqc_%J.out" \
        --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} \
        "${moduledir}/2C-fastqc-trim.sh"
) &
prep_pids+=($!)

run 3A-star_index_genome \
    --error="${log}/3A-star_index_%J.err" --output="${log}/3A-star_index_%J.out" \
    --job-name=${NAME} --partition=${PART} "${moduledir}/3A-star_index_genome.sh" &
prep_pids+=($!)

# Gate: mapping needs trimmed reads (2B) and the index (3A)
prep_failed=0
for pid in "${prep_pids[@]}"; do
    wait "$pid" || prep_failed=1
done
[[ "$prep_failed" -ne 0 ]] && { echo "prep failed - not starting mapping" >&2; exit 1; }

###### Step 3B: star mapping
run 3B-star_mapping_array \
    --error="${log}/3B-star_map_%J.err" --output="${log}/3B-star_map_%J.out" \
    --array="0-${sample_number}%14" --job-name=${NAME} --partition=${PART} \
    "${moduledir}/3B-star_mapping_array.sh"

###### Step 4: MapQC (produces the markdup/rmdup BAMs used by step 5)
run 4-mapqc_array \
    --error="${log}/4-mapqc_%J.err" --output="${log}/4-mapqc_%J.out" \
    --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} \
    "${moduledir}/4-mapqc_array.sh"

###### Step 5: featurecount
run 5-featurecount_array \
    --error="${log}/5-featurecount_%J.err" --output="${log}/5-featurecount_%J.out" \
    --array="0-${sample_number}%20" --job-name=${NAME} --partition=${PART} \
    "${moduledir}/5-featurecount_array.sh"

###### Step 6: SARTools
run 6-SARTools \
    --error="${log}/6-sartools_%J.err" --output="${log}/6-sartools_%J.out" \
    --job-name=${NAME} --partition=${PART} "${moduledir}/6-SARTools.sh"

###### Step X: MultiQC report
run X-multiqc \
    --error="${log}/multiqc_%J.err" --output="${log}/multiqc_%J.out" \
    --job-name=${NAME} --partition=${PART} "${moduledir}/X-multiqc.sh"

echo "$(date '+%F %T') PIPELINE COMPLETE"
