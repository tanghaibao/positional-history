#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Batch process the blast filter and calculate synteny score for a group of species
make it ready for synteny_graph
"""

import sys
from subprocess import Popen
from string import Template

def sh(cmd, go=True):
    print >>sys.stderr, cmd
    if go:
        proc = Popen(cmd, shell=True)
        out, err = proc.communicate()

# paths for data
S="athaliana lyrata papaya peach grape".split()
BDIR="/home/bao/blast/blastz"
DDIR="bed"
QDIR="/home/bao/projects/quota-alignment"
RDIR="results"


# run the batch commands
len_S = len(S)
for species in S[1:]:
    a, b = "athaliana", species
    if a > b: a, b = b, a
    cmd1 = "blast_to_raw.py ${BDIR}/${a}_${b}.blastz --write-filtered-blast --tandem_Nmax 10 --qbed=${DDIR}/$a.bed --sbed=${DDIR}/$b.bed >${a}_${b}.blastz.filtered"
    cmd2 = "synteny_score.py ${a}_${b}.blastz.filtered --qbed=${DDIR}/${a}.bed --sbed=${DDIR}/$b.nolocaldups.bed >${RDIR}/${a}_${b}.synteny_score"
    cmd1 = Template(cmd1).substitute(locals())
    cmd2 = Template(cmd2).substitute(locals())
    go = True 
    #sh(cmd1, go=go)
    sh(cmd2, go=go)

