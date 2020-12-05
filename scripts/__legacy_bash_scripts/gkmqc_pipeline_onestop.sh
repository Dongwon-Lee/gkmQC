#!/bin/bash
#SBATCH --job-name=gqc_1stop
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=bch-compute
#SBATCH --mem=4G
##SBATCH --mail-type=end
##SBATCH --mail-user=seongkyu.han@childrens.harvard.edu
module load anaconda3

PROJ_DIR=/home/ch218391/PROJECTS/gkmqc
SCRIPT_DIR=${PROJ_DIR}/scripts
BED_DIR=$PWD # Generalize for all cases
#BED_DIR=${PROJ_DIR}/analysis/encode_dnase_seq_bed

INPUTF=$1 # narrowPeak file
GENOME=$2 # hg38
PREFIX=${INPUTF%_peaks.narrowPeak}

# Make/goto and workdir
PREPROC_DIR=${BED_DIR}/${PREFIX}.gkmqc
mkdir ${PREPROC_DIR}
cd ${PREPROC_DIR}
cp ../${INPUTF} .

# 1. Preprocessing scripts
echo "first stage: generating spliited data..."
PREPROC_SCRIPT=${SCRIPT_DIR}/gkmqc_pipeline_preproc.sh
sh ${PREPROC_SCRIPT} ${INPUTF} ${2}

cd ${PREPROC_DIR}
# 2. Generating splitted data
echo "second stage: 5-fold cross-validation"
CV_SCRIPT=${SCRIPT_DIR}/gkmqc_pipeline_5cv.sh
sh ${CV_SCRIPT} ${INPUTF}
