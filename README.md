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
  -h, --help          Show this help message
  -h, --help          Show this help message
```
 **Note:**
- You can run the pipeline multiple times simultaneously with different raw reads, simply repeat the installation process in a different directory and `./deploy` with a different run names identifier name.
- You can manually reconfigure slurm parameters as per your HPC system (e.g memory, CPUs) by going through indivudal scripts in `modules` directory.
- All the relevent outputs will be stored in `outdir` folder, and outputs for every individual steps in the pipeline can be found in `workdir`.

Prof Peter Kille - kille@cardiff.cf.ac.uk
