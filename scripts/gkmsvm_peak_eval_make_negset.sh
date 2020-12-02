#!/bin/bash
#SBATCH --job-name=mknegset
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=bch-compute
#SBATCH --mem=4G
# #SBATCH --mail-type=end
# #SBATCH --mail-user=dongwon.lee@childrens.harvard.edu

set -o errexit
set -o nounset

POSF=$1 # positive bed file (example: hrtaa_atac_r1.e300.qc.top1.bed)
GENOME=$2 # genome version

SCRIPTDIR=/home/ch218391/PROJECTS/gkmqc/scripts
NULLSEQDIR=/lab-share/Neph-Sampson-e2/Public/dlee/nullseq_indice/$GENOME
GENOMESEQDIR=/lab-share/Neph-Sampson-e2/Public/dlee/genomes/$GENOME

PREFIX=${POSF%.bed}
NEGF=${PREFIX}.nr1.bed
NEGF_FASTA=${PREFIX}.nr1.fa

# 1. generate a negative set
if [[ -e $NEGF ]]; then
    NPB=`cat $POSF |wc -l`
    NNB=`cat $NEGF |wc -l`
    if [[ $NPB -eq $NNB ]]; then
        echo "skip making $NEGF"
    else
        python ${SCRIPTDIR}/nullseq_generate.py -o $NEGF $POSF $GENOME $NULLSEQDIR
    fi
else
    python ${SCRIPTDIR}/nullseq_generate.py -o $NEGF $POSF $GENOME $NULLSEQDIR
fi

# 2. make a fasta file
if [[ -e $NEGF_FASTA ]]; then
    N=`cat $NEGF| wc -l`
    N2=$(( $N * 2 ))
    S2=`cat $NEGF_FASTA|wc -l`
    if [[ $N2 -eq $S2 ]]; then
        echo "skip making $NEGF_FASTA"
    else
        python ${SCRIPTDIR}/fetchseqs.py -d $GENOMESEQDIR $NEGF $NEGF_FASTA
    fi
else
    python ${SCRIPTDIR}/fetchseqs.py -d $GENOMESEQDIR $NEGF $NEGF_FASTA
fi
