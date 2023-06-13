#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J launch_cromwell
#SBATCH --mem=64000
#SBATCH -p day
#SBATCH -t 1400


java -Dconfig.file=slurm.conf -jar \
/home/ml2529/palmer_scratch/SMuRF/cromwell-56.jar run \
gargamel.wdl \
-i fkrp.json \
-o cromwell.options

