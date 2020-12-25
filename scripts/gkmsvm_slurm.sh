#!/bin/bash
#SBATCH --job-name=gkmsvm_wrapper
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=bch-compute
#SBATCH --mem=10G
#SBATCH --export=NONE

#set -o errexit
#set -o nounset
source activate gkmqc
$@