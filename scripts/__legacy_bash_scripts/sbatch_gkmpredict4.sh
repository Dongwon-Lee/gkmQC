#!/bin/bash
#Wrapper script for sbatch with 4-threads

#SBATCH --job-name=gkmpred4
#SBATCH --time=11:59:59
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu_short
#SBATCH --mem=20G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dwlee@jhu.edu

set -o errexit
set -o nounset

~/program/lsgkm/bin/gkmpredict -T 4 "$@"
