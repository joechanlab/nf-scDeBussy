#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --job-name=nextflow
#SBATCH -o log/%A.out
#SBATCH -e log/%A.err
#SBATCH --time=06:00:00

module load miniforge3
source activate /usersoftware/chanj3/nextflow
export PYTHONPATH=$PYTHONPATH:/home/wangm10/:/usersoftware/chanj3/tslearn
nextflow run ./main.nf -resume -profile iris -params-file ./data/params.yml -w "/scratch/chanj3/wangm10/scDeBussy_work"
