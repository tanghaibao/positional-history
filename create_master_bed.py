#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from bed_utils import Bed
from orgs import species

if __name__ == '__main__':
    fw = open("bed/all.bed", "w")
    for s in species:
        bed = Bed("bed/%s.bed" % s)
        for b in bed:
            seqid = b.seqid
            if s=="athaliana": 
                seqid = seqid.replace("Chr", "")
            b.seqid = "%s:%s" % (s, seqid)
            print >>fw, b
