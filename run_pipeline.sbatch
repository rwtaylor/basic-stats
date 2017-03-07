#!/bin/bash
#SBATCH --job-name=bs-170208-01
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=dpetrov
#################

./nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace -with-dag dag.png
