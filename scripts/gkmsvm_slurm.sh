#!/bin/bash
#SBATCH --job-name=gkmsvm_wrapper
#SBATCH --time=23:59:59
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=bch-compute
#SBATCH --mem=15G

#set -o errexit
#set -o nounset
source ~/program/miniforge3/bin/activate scatac
$@