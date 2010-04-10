#!/bin/bash

set -v

S="lyrata papaya poplar medicago grape"
BDIR=/home/bao/blast/results
DDIR=bed
QDIR=/home/bao/projects/quota-alignment
RDIR=results

i=peach
${QDIR}/scripts/blast_to_raw.py ${BDIR}/athaliana_${i}.blastp --tandem_Nmax 10 --cscore .7 --filter_repeats  --qbed=${DDIR}/athaliana.bed --sbed=${DDIR}/${i}.bed >athaliana_${i}.blastp.filtered
${QDIR}/scripts/synteny_score.py athaliana_${i}.blastp.filtered --qbed=${DDIR}/athaliana.bed --sbed=${DDIR}/${i}.nolocaldups.bed >${RDIR}/athaliana_${i}.synteny_score 

#for i in ${S}; do
#    ${QDIR}/scripts/synteny_score.py athaliana_${i}.blastp.filtered --qbed=${DDIR}/athaliana.bed --sbed=${DDIR}/${i}.nolocaldups.bed >${RDIR}/athaliana_${i}.synteny_score 
#done
