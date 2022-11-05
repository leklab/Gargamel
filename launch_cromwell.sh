#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J launch_cromwell
#SBATCH --mem=16000

java -Dconfig.file=slurm.conf -jar \
/gpfs/ycga/project/lek/shared/tools/jars/cromwell-56.jar run \
smurf.wdl \
-i smurf.json \
-o cromwell.options

