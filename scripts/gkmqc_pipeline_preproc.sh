#!/bin/bash
#SBATCH --job-name=gkmqc
#SBATCH --time=8:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=bch-compute
#SBATCH --mem=4G
##SBATCH --mail-type=end
##SBATCH --mail-user=seongkyu.han@childrens.harvard.edu

INPUTF=$1 # MACS narrowPeak file (assuming that it ends with _peaks.narrowPeak).
GENOME=$2 # hg38/mm10/etc...

EXTLEN=300
NPEAKS=5000
SCRIPTDIR=/home/ch218391/PROJECTS/gkmqc/scripts

PREFIX=${INPUTF%_peaks.narrowPeak}.e${EXTLEN}

POSF=${PREFIX}.qc.bed

MAKE_QC_POSSET_SCRIPT=${SCRIPTDIR}/gkmsvm_peak_eval_make_qc_posset.sh
SPLIT_POSSET_SCRIPT=${SCRIPTDIR}/gkmsvm_peak_eval_split_posset.sh
MAKE_NEGSET_SCRIPT=${SCRIPTDIR}/gkmsvm_peak_eval_make_negset.sh

# 1. QC and make a positive set
sh ${MAKE_QC_POSSET_SCRIPT} ${INPUTF} ${EXTLEN} ${GENOME}

# 2. split the positive set by p-value 
sh ${SPLIT_POSSET_SCRIPT} ${POSF} ${NPEAKS}

# 3. generate a negative set for each of the splitted positive sets (1..19)
for fn in ${PREFIX}.qc.top?.bed ${PREFIX}.qc.top1?.bed; do
    if [[ -e $fn ]]; then
        #3.1. get a positive set
        python ${SCRIPTDIR}/select_seqs.py ${PREFIX}.qc.fa ${fn} >${fn%.bed}.fa

        #3.2. get a negative set
        sh ${MAKE_NEGSET_SCRIPT} ${fn} ${GENOME}
    fi
done
