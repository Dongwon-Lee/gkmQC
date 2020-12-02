#!/bin/bash
#SBATCH --job-name=mkposqc
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=bch-compute
#SBATCH --mem=4G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dongwon.lee@childrens.harvard.edu

set -o errexit
set -o nounset

INPUTF=$1 # MACS narrowPeak file (assuming that it ends with _peaks.narrowPeak).
EXTLEN=$2 # PEAKLEN/2 (ie. 150 for 300bp, 300 for 600bp)
GENOME=$3 # genome version
PREFIX=${INPUTF%_peaks.narrowPeak}.e${EXTLEN}

SCRIPTDIR=/home/ch218391/PROJECTS/gkmqc/scripts
NULLSEQDIR=/lab-share/Neph-Sampson-e2/Public/dlee/nullseq_indice/$GENOME
GENOMESEQDIR=/lab-share/Neph-Sampson-e2/Public/dlee/genomes/$GENOME

POSF0=${PREFIX}.bed
POSF0_PROF=${PREFIX}.prof

POSF=${PREFIX}.qc.bed
NEGF=${PREFIX}.qc.nr1.bed

POSF_FASTA=${PREFIX}.qc.fa
NEGF_FASTA=${PREFIX}.qc.nr1.fa

# 1. make fixed length peaks
if [[ -e $POSF0 ]]; then
    echo "skip making $POSF0"
else
    awk -v OFS="\t" -v SHFT=$EXTLEN \
    '$1 != "chrM" && $1 != "chrY" && $1 !~ /_random/ && $1 !~ /chrUn/ && $1 !~ /hap/ && $1 !~ /_alt/ && $2+$10-SHFT > 0 {
    summit=$2+$10;
    print $1,summit-SHFT,summit+SHFT,$4,$7}' $INPUTF >$POSF0 ## CHANGED FOR HOTSPOT2 INPUT ## $4,$8 -> mac2, $4,$7 -> hotspot2
fi

# 2. calculate profiles of the fixed length peaks
NB=`cat $POSF0|wc -l`
if [[ -e $POSF0_PROF ]]; then
    NP=`cat $POSF0_PROF|wc -l`
    if [[ $NB -eq $NP ]]; then
        echo "skip making $POSF0_PROF"
    else
        python ${SCRIPTDIR}/make_profile2.py $POSF0 $GENOME $NULLSEQDIR $POSF0_PROF
    fi
else
    python ${SCRIPTDIR}/make_profile2.py $POSF0 $GENOME $NULLSEQDIR $POSF0_PROF
fi

# 3. remove peaks with >1% of N bases & >70% of repeats
if [[ -e $POSF ]]; then
    echo "skip making $POSF"
else
    paste $POSF0_PROF $POSF0 | awk '$4<=0.7 && $5<=0.01' |cut -f 6- >$POSF
fi

##########################################################

# 5. make fasta files
if [[ -e $POSF_FASTA ]]; then
    N=`cat $POSF| wc -l`
    N2=$(( $N * 2 ))
    S2=`cat $POSF_FASTA|wc -l`
    if [[ $N2 -eq $S2 ]]; then
        echo "skip making $POSF_FASTA"
    else
        python ${SCRIPTDIR}/fetchseqs.py -d $GENOMESEQDIR $POSF $POSF_FASTA
    fi
else
    python ${SCRIPTDIR}/fetchseqs.py -d $GENOMESEQDIR $POSF $POSF_FASTA
fi
