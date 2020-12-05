#!/bin/bash
#Wrapper script for sbatch with 4-threads

#SBATCH --job-name=gkmtrain4
#SBATCH --time=11:59:59
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=bch-compute
#SBATCH --mem=15G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dongwon.lee@childrens.harvard.edu

set -o errexit
set -o nounset

~/program/lsgkm/bin/gkmtrain -T 4 -m 10000 "$@"
