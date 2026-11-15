# RNASeq analysis from Reads to DEGs
A SLURM pipeline designed for illumina data using (star)[https://github.com/alexdobin/STAR/tree/master] mapping algorithum, (subread)[https://subread.sourceforge.net/] and to generate feature counts and [SARTools](https://github.com/PF2-pasteur-fr/SARTools) for DEG geneartion. The processes includes quality assessment, pre-processing, kraken, bracken, metaphylan with alpha diversity metrics being calculated and krona plots being generated. This pipeline is optimised explicitly for deployment on high-performance computing (HPC) clusters and executed _via_ Slurm Workload Manager.

## Key Features

- **Read Quality Assessment**: Ensures high-quality data processing with integrated quality checks.
- **Library indexing and Mapping**: Star.
- **Indentification of duplicates markup/removal files **: Picard.
- **FeatureCount generation**: Subread - note by default this is configured for eukaryotes and exon based annotations.
- **Generation of DEGs**: Generates DEGs using SARTools pipeline and saving of DESeq2 object


## Installation

1. Install the metagenome_slurm resources into your HPC cluster directory in which you will be performing the assembly:  

```
git clone https://github.com/Peter-Kille/transcript_reads2DEGs_slurm.git
```

2. Put the raw reads and metadata in `source_data` folder and genome (.fasta) and annotation (.gtf) in `genome` folder.  

3. Run the pipeline using `./deploy.sh -n [NAME] -p [PARTITION] -m [METADATA] -t [TREATMENT] -r [REFERENCE]`  

4. you can you './deploy.sh -h' for help (see below)

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
## How the run is controlled (v2: driver job, no singleton chain)

`./deploy.sh` no longer submits every stage up front chained with `sbatch -d singleton`. It validates your options, then submits itself as a small **driver** job (1 CPU, 1 GB, default 7 days on `epyc_long`) and exits. The driver submits one stage at a time with `sbatch --wait`, so only the driver and the current stage are ever in the queue (this keeps array tasks within `MaxSubmitJobs`). The design follows the [metagenomic_slurm pipeline guide](https://github.com/Peter-Kille/metagenomic_slurm/blob/main/PIPELINE_GUIDE.md).

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
- The previous singleton-based version is preserved as release **v1**.
 **Note:**
- You can run the pipeline multiple times simultaneously with different raw reads, simply repeat the installation process in a different directory and `./deploy` with a different run names identifier name.
- You can manually reconfigure slurm parameters as per your HPC system (e.g memory, CPUs) by going through indivudal scripts in `modules` directory.
- All the relevent outputs will be stored in `outdir` folder, and outputs for every individual steps in the pipeline can be found in `workdir`.

Prof Peter Kille - kille@cardiff.cf.ac.uk
