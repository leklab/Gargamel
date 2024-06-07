**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies. The manuscript is currently available on [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3). This workflow involves a 2-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated with short-read Next-generation sequencing. 

We have organized scripts used for SMuRF across three repositories:    
* [**Balthazar repository**](https://github.com/leklab/Balthazar) is used to create the oligo library needed for the saturation mutagenesis as well as perform QC for the resulting construct.    
* **Gagarmel repository (this repository)** is used to process the raw data from NGS to quantify the enrichment of the variants.    
* [**Azrael repository**](https://github.com/leklab/Azrael) is used to generate and analyze the functional scores, and plot the results.    

## Requirements

* Java 1.8
* Python (tested on v3.8.6)
* Cromwell (tested on v56)
* samtools (tested on v1.16)
* bwa (tested on v0.7.17)

### Python libraries required
```
pysam (test on v0.16.0.1)
```

## Step 1: Loading required libraries
Tested on McCleary cluster at Yale
```
module load Java/1.8.345
module load Python/3.8.6-GCCcore-10.2.0
module load SAMtools/1.16-GCCcore-10.2.0
module load BWA/0.7.17-GCC-10.2.0
module load Pysam/0.16.0.1-GCC-10.2.0
```

## Step 2: Preparation of input files

* Change the output directory paths in the **cromwell.options** file
* Change the paths in **large1.json** and **fkrp.json** to reference the absolute paths for the scripts and reference_files folders.
* Change the input files to reference the demultiplexed FASTQ files and the block demarcations

## Step 3: Submit the job to the Slurm job manager
```
sbatch launch_cromwell.sh
```
Output files will be located in the directory defined in the **cromwell.options** file
