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
from orgs import species as S
BDIR="blasts"
DDIR="bed"
QDIR="/home/bao/projects/quota-alignment"
RDIR="results"

BLAST_NB="/usr/local/src/bio-pipeline/blast_nearby/blast_nearby.py"
BLAST_NB_CMD="bl2seq -i %(query_fasta)s -j %(subject_fasta)s -p blastn -e 0.05 -W 15 -D 1 "
FDIR='fastas/'

cmds = ["./blast_to_raw.py ${BDIR}/${a}_${b}.blastz --write-filtered-blast --tandem_Nmax 10 --qbed=${DDIR}/$a.bed --sbed=${DDIR}/$b.bed > ${BDIR}/${a}_${b}.blastz.filtered",
        "./synteny_score.py ${BDIR}/${a}_${b}.blastz.filtered --qbed=${DDIR}/${a}.bed --sbed=${DDIR}/${b}.nolocaldups.bed >${RDIR}/${a}_${b}.synteny_score",
        "./extract_gray_anchors.py --qbed=${DDIR}/${a}.bed --qfasta ${FDIR}/${a}.fa  --sbed=${DDIR}/${b}.nolocaldups.bed --sfasta ${FDIR}/${b}.fa --pad 15000 ${RDIR}/${a}_${b}.synteny_score --out-prefix ${RDIR}/${a}_${b}",
        "${BLAST_NB} --anchors ${RDIR}/${a}_${b}.anchors --dist 0 ${FDIR}/${a}.fa ${FDIR}/${b}.fa --cmd '${BLAST_NB_CMD}'  > ${BDIR}/${a}_${b}.nearby.blast",
        "./count_hits_near_flankers.py --missing ${RDIR}/${a}_${b}.missing --qbed bed/${a}.bed --sbed bed/${b}.nolocaldups.bed --dist 15000 --blast ${BDIR}/${a}_${b}.nearby.blast ${RDIR}/${a}_${b}.synteny_score > ${RDIR}/${a}_${b}.nearby.synteny_score"
       ]

# run the batch commands
len_S = len(S)

# modify range of species to run certain species
for species in ["castor"]:
    a, b = "athaliana", species
    if a > b: a, b = b, a
    print "\n" + ("=" * len(b))
    print b
    print "=" * len(b)
    go = True
    # modify range of cms to run part of pipeline
    for cmd in cmds:
        cmd = Template(cmd).substitute(locals())
        sh(cmd, go=go)

# final command to execute (yay)
cmd = "python assemble_synteny.py"
sh(cmd)
