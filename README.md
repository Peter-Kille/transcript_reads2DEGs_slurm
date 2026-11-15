# RNASeq analysis from Reads to DEGs
A SLURM pipeline designed for Illumina paired-end RNA-seq data. It uses [STAR](https://github.com/alexdobin/STAR/tree/master) for mapping, [Subread featureCounts](https://subread.sourceforge.net/) to generate feature counts, and [SARTools](https://github.com/PF2-pasteur-fr/SARTools) (DESeq2) for differential expression (DEG) generation. The process includes quality assessment (FastQC, fastp), mapping QC (Picard, RustQC), feature counting, DEG analysis and a MultiQC report. The pipeline is optimised for deployment on high-performance computing (HPC) clusters and is executed _via_ the Slurm Workload Manager.

## Versions

| Version | What it is | Get it |
|---|---|---|
| **main** (current) | Driver-job design: one stage submitted at a time with `sbatch --wait`, done markers, resumable (see [How the run is controlled](#how-the-run-is-controlled-current-version)). | `git clone https://github.com/Peter-Kille/transcript_reads2DEGs_slurm.git` |
| **v1** | The original pipeline. See [Version 1](#version-1-original-singleton-pipeline) below. | `git clone --branch v1 https://github.com/Peter-Kille/transcript_reads2DEGs_slurm.git` |

## Key Features

- **Read Quality Assessment**: FastQC on raw and trimmed reads, fastp trimming.
- **Library indexing and Mapping**: STAR.
- **Identification of duplicates markup/removal files**: Picard.
- **FeatureCount generation**: Subread - note by default this is configured for eukaryotes and exon based annotations.
- **Generation of DEGs**: Generates DEGs using the SARTools pipeline and saves the DESeq2 object.

## Installation

1. Clone the pipeline into the HPC cluster directory in which you will run the analysis (current version shown; for v1 see [below](#how-to-clone-version-1)):

```
git clone https://github.com/Peter-Kille/transcript_reads2DEGs_slurm.git
cd transcript_reads2DEGs_slurm
```

2. Put the raw reads and metadata in the `source_data` folder and the genome (.fasta) and annotation (.gtf) in the `genome` folder.

3. Run the pipeline using `./deploy.sh -n [NAME] -p [PARTITION] -m [METADATA] -t [TREATMENT] -r [REFERENCE]`

4. You can use `./deploy.sh -h` for help (see below).

## Available displayed arguments:
```
./deploy.sh -h

Usage: ./deploy.sh -n [NAME] -p [PARTITION] -m [METADATA] -t [TREATMENT] -r [REFERENCE]

Copy your paired end sequence files as .fastq.gz files in the source_data directory. Also
copy a metadata file into the source_data directory.  An example tab deliminated (.tsv)
metadata file is provided.  First column should have sample name, column 2 should be
full name forward read and third column the reverse read. Column 4 should be your treatment.
The treatment column should contain refernce term - avoid all symbols or spaces
when creating the descriptions of your sample treatments.

  REQUIRED: place an uncompressed genome file into a genome directory ensuring it has suffix .fasta
  REQUIRED: place an uncompressed annotation file into a genome directory ensuring it has suffix .gtf

Options:
  -n, --name          REQUIRED: Run name or deployment name - should be unique
  -p, --partition     REQUIRED: Avalible partition / hpc queue (epyc, defq, jumbo, epyc_ssd)
  -m, --metadata      REQUIRED: metadata file, tsv col=sample name, col2=read1, col3=read2
  -t, --treatment     REQUIRED: treatment for which to derive DEGs
  -r, --ref           REQUIRED: reference condition
  -w, --work          Optional: working dir - default is current dir /work/
  -D, --driver-partition  Optional: partition for the driver job (default epyc_long)
  -T, --driver-time   Optional: driver walltime (default 7-00:00:00)
  -F, --foreground    Optional: run the stages in this shell instead of a driver job (use in tmux)
  -h, --help          Show this help message
```
The `-D`, `-T` and `-F` options exist only in the current version; v1 does not have them.

## How the run is controlled (current version)

`./deploy.sh` submits itself as a small **driver** job (1 CPU, 1 GB, default 7 days on `epyc_long`) and exits. The driver submits one stage at a time with `sbatch --wait`, so only the driver and the current stage are ever in the queue (this keeps array tasks within `MaxSubmitJobs`). The design follows the [metagenomic_slurm pipeline guide](https://github.com/Peter-Kille/metagenomic_slurm/blob/main/PIPELINE_GUIDE.md).

```
1-preprocess
   -> [ 2A fastqc -> 2B fastp -> 2C trim fastqc ]  in parallel with  [ 3A STAR index ]
   -> wait for both -> 3B STAR mapping -> 4 mapqc -> 5 featureCounts -> 6 SARTools -> X MultiQC
```

- **Done markers:** each finished stage writes `work/.done_<stage>`. If a stage fails the driver stops; fix the problem and run the *same* command again - finished stages print `SKIP`. To force a stage to rerun delete its marker and the markers of every later stage.
- **Follow progress:** `tail -f work/0_log/0_driver_<jobid>.out` (lines read SUBMIT / DONE / SKIP / FAILED); stage logs are `work/0_log/<stage>_<jobid>.out/.err`.
- **Cancel:** `scancel <driver_id>`, then `scancel -u $USER -n <NAME>` for the stage in flight (cancelling the driver does not stop the running stage).
- **Command history:** every launch is appended to `work/0_log/0_execution_commandline.txt`.
- **Debugging / tmux:** use `-F` to skip the driver and run the stages in the current shell.
- Use a new `-n` name (and a clean `work/`) for a fresh run, otherwise done markers will skip stages you meant to repeat.

## Version 1 (original singleton pipeline)

**What it is.** Release [v1](https://github.com/Peter-Kille/transcript_reads2DEGs_slurm/releases/tag/v1) is the pipeline as it stood before the driver-job upgrade (commit `86082ce`, "star array job limit 14"). It has the same stages and the same analysis steps as the current version; only the way the stages are submitted differs:

- `deploy.sh` submits **every stage at once** from the login node, each chained with `sbatch -d singleton`. Slurm's singleton dependency runs jobs sharing the same job name (`-n`) and user one after another, so the stages run strictly in series: no stage overlaps, including the QC steps and STAR index build.
- There is no driver job, no `sbatch --wait`, no done markers and no `-D`, `-T` or `-F` options.
- Because all stages (and every array task) are queued upfront, large sample sets can hit the cluster's `MaxSubmitJobs` limit, and after a failure there is no automatic resume: you must clean up `work/` and resubmit by hand.
- Known issues fixed in the current version: the `--metadata` long option was misspelt (`--metdata`), a single-sample metadata file was rejected as "no files", the step 5 array tasks all copied the shared output folder at once, and MultiQC failed on a rerun when `outdir/multiqc` already existed.

Use v1 if you need to reproduce results exactly as they were generated with it. Otherwise use the current version.

### How to clone version 1

Clone the tagged release straight into a new folder (HTTPS):

```
git clone --branch v1 https://github.com/Peter-Kille/transcript_reads2DEGs_slurm.git transcript_reads2DEGs_slurm_v1
cd transcript_reads2DEGs_slurm_v1
```

If you use an SSH key with GitHub:

```
git clone --branch v1 git@github.com:Peter-Kille/transcript_reads2DEGs_slurm.git transcript_reads2DEGs_slurm_v1
cd transcript_reads2DEGs_slurm_v1
```

Git will say you are in a "detached HEAD" state. That is expected for a tag and is fine for running the pipeline; do not commit changes there without creating a branch (`git switch -c my-branch`).

Already cloned the current version? Switch the existing clone to v1 and back:

```
git fetch --tags
git checkout v1          # original singleton version
git checkout main        # back to the current version
```

Check which version you have:

```
git describe --tags      # prints v1 when on the release
git log -1 --oneline     # v1 is 86082ce "star array job limit 14"
```

No git? Download the v1 source as an archive from the [release page](https://github.com/Peter-Kille/transcript_reads2DEGs_slurm/releases/tag/v1) (Source code, tar.gz or zip), or with:

```
wget -O transcript_reads2DEGs_slurm_v1.tar.gz https://github.com/Peter-Kille/transcript_reads2DEGs_slurm/archive/refs/tags/v1.tar.gz
mkdir transcript_reads2DEGs_slurm_v1 && tar -xzf transcript_reads2DEGs_slurm_v1.tar.gz -C transcript_reads2DEGs_slurm_v1 --strip-components=1
```

### How to run version 1

Setup is the same as the current version (reads and metadata in `source_data`, genome `.fasta` and `.gtf` in `genome`). Then:

```
./deploy.sh -n [NAME] -p [PARTITION] -m [METADATA] -t [TREATMENT] -r [REFERENCE]
```

This queues all stages and returns; follow progress with `squeue -u $USER -n [NAME]` and the logs in `work/0_log/`. Use a unique `-n` for each run, since the singleton dependency is keyed on the job name.

## Notes
- You can run the pipeline multiple times simultaneously with different raw reads: repeat the installation in a different directory and run `./deploy.sh` with a different run name.
- You can manually reconfigure Slurm parameters as per your HPC system (e.g. memory, CPUs) by editing the individual scripts in the `modules` directory.
- All the relevant outputs are stored in the `outdir` folder, and outputs for every individual step can be found in `work`.

Prof Peter Kille - kille@cardiff.cf.ac.uk
