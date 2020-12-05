#!/bin/bash

INPUTF=$1 # MACS narrowPeak file (assuming that it ends with _peaks.narrowPeak).

EXTLEN=300
SCRIPTDIR=/home/ch218391/PROJECTS/gkmqc/scripts

PREFIX=${INPUTF%_peaks.narrowPeak}.e${EXTLEN}
GKMSVM_SCRIPT=${SCRIPTDIR}/sbatch_gkmtrain4.sh
T=4
L=10
K=6
D=3
NCV=5

# 3. train with 5-CV for each of the splitted positive sets... up to 19 (1..19)
for fn in ${PREFIX}.qc.top?.fa ${PREFIX}.qc.top1?.fa; do
    if [[ -e $fn ]]; then
        echo $fn
        # positive set
        POSF=${fn}
        NEGF=${POSF%.fa}.nr1.fa

        sbatch ${GKMSVM_SCRIPT} -t ${T} -l ${L} -k ${K} -d ${D} -x ${NCV} ${POSF} ${NEGF} gkmsvm.4.10.6.3.${POSF%.fa}
        sleep 0.1
    fi
done
