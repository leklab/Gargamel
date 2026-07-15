# Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)

**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies.

The original SMuRF study, **“Saturation mutagenesis-reinforced functional assays for disease-related genes,”** is published in [Cell](https://www.cell.com/cell/fulltext/S0092-8674%2824%2900976-0).

This workflow uses two-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated using short-read next-generation sequencing.

## Repositories

Scripts used for the original SMuRF workflow are organized across three repositories:

* The [**Balthazar repository**](https://github.com/leklab/Balthazar) is used to create the oligonucleotide library required for saturation mutagenesis and to perform quality control on the resulting construct.
* The **Gargamel repository (this repository)** is used to process raw next-generation sequencing data and quantify variant enrichment.
* The [**Azrael repository**](https://github.com/leklab/Azrael) is used to generate and analyze functional scores and plot the results.

An additional sister study applies the SMuRF workflow to the **SGCA** gene:

* The [**Denisa repository**](https://github.com/leklab/denisa) contains all scripts and resources relevant to the SGCA study, including data processing, membrane localization score generation, confidence classification, variant annotation, and downstream analyses.

## Requirements

* Java 1.8
* Python, tested on v3.8.6
* Cromwell, tested on v56
* SAMtools, tested on v1.16
* BWA, tested on v0.7.17

### Required Python libraries

```text
pysam, tested on v0.16.0.1
```

## Step 1: Load the required software

The following module versions were tested on the McCleary cluster at Yale:

```bash
module load Java/1.8.345
module load Python/3.8.6-GCCcore-10.2.0
module load SAMtools/1.16-GCCcore-10.2.0
module load BWA/0.7.17-GCC-10.2.0
module load Pysam/0.16.0.1-GCC-10.2.0
```

## Step 2: Prepare the input files

* Change the output directory paths in the `cromwell.options` file.
* Change the paths in `large1.json` and `fkrp.json` so that they reference the absolute paths of the `scripts` and `reference_files` directories.
* Change the input-file paths so that they reference the demultiplexed FASTQ files.
* Specify the appropriate block demarcations for the gene being processed.

## Step 3: Submit the job to the Slurm job manager

```bash
sbatch launch_cromwell.sh
```

Output files will be written to the directory specified in the `cromwell.options` file.
