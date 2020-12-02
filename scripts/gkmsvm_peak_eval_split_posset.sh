#!/bin/bash

#SBATCH --job-name=splitpos
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=bch-compute
#SBATCH --mem=4G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dongwon.lee@childrens.harvard.edu

set -o errexit
set -o nounset

BEDF=$1 # bed file after QC (example: hrtaa_atac_r1.e300.qc.bed)
npeaks=$2 # number of peaks per test: 10000
#XXX: assume that number of peaks should be more than $npeaks

ntot=`cat $BEDF |wc -l`
ntests=$(( ($ntot + $npeaks/2) / $npeaks ))

sort -gr -k5,5 $BEDF >${BEDF}.tmp.sorted

for i in `seq $ntests`; do
    SKIPN=$(( ($i - 1) * $npeaks + 1 ))
    if [[ $i -eq $ntests ]]; then # handling the last case. No need to take $npeaks
        tail -n +$SKIPN ${BEDF}.tmp.sorted |sortBed >${BEDF%.bed}.top${i}.bed
    else
        tail -n +$SKIPN ${BEDF}.tmp.sorted |head -n $npeaks |sortBed >${BEDF%.bed}.top${i}.bed
    fi
done

rm -f ${BEDF}.tmp.sorted
