#!/bin/bash
#Wrapper script for sbatch with 4-threads

#SBATCH --job-name=gkmweights4
#SBATCH --time=2:00:00
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
~/program/lsgkm/bin/gkmpredict ~/PROJECTS/gkmqc/data/10mers_nr.fa ${11}.model.txt ${11}.10mers.txt
