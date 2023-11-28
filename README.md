**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies. The manuscript is currently available on BioRxiv: https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3. 


This workflow involves a 2-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated with short-read Next-generation sequencing. We have organized scripts used for SMuRF across three repositories:    

**Balthazar repository** is used to create the oligo library needed for the saturation mutagenesis as well as perform QC for the resulting construct.    
**Gagarmel repository (this repository)** is used to process the raw data from NGS to quantify the enrichment of the variants.    
**Azrael repository** is used to generate and analyze the functional scores, and plot the results.    


## Requirements

* Java 1.8
* Python (tested on v3.8.6)
* Cromwell (tested on v56)
* samtools (tested on v1.16)
* bwa (tested on v0.7.17)

### Python libraries required
```
pysam (test ov v0.16.0.1)
```
