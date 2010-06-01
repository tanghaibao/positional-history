#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import csv
import itertools
import collections
from bed_utils import Bed

"""
import bitly
api = bitly.Api(login="tanghaibao", apikey="R_51b7c08fbb167cc934ca6ae7fd9e8b2d")
def shorten(a, url):
    return a.shorten(url)
"""
species = "athaliana lyrata grape papaya peach".split()
DSIDS = [42029, 39219, 43318, 563, 42478]
DSIDS = dict(zip(species, DSIDS))


class SyntenyLine(object):
    # GSVIVT01012261001       POPTR_0008s09790        G       5       80000   -
    def __init__(self, line):
        args = line.strip().split("\t")
        self.query = args[0]
        self.anchor = args[1]
        self.gray = args[2]
        self.score = int(args[3])
        self.dr = int(args[4])
        self.rev = args[5]
    def __str__(self):
        return "\t".join(str(x) for x in (self.query, self.anchor, self.gray, self.score, self.dr, self.rev))
    __repr__=__str__

quota = dict(grape=1, papaya=1, poplar=2, medicago=2, lyrata=1, peach=1)


def attach_species(rec, s):
    fp = file("results/athaliana_%s.synteny_score" % s)
    q = quota[s] # get top n hits (ranked by synteny score)
    qrec = collections.defaultdict(list)
    for row in fp:
        atoms = row.split()
        if atoms[1]=="na": continue
        s = SyntenyLine(row)
        query = s.query
        if not query.startswith("AT"): continue
        if len(qrec[query]) < q:
            qrec[query].append(s)

    for q, v in qrec.iteritems():
        rec[q].append(sorted(v, key=lambda x: x.anchor))


def build_gevo_link(query, group, i, bed, coge_pos):

    global species
    # make query region, +-20 genes
    seqid, start = bed[i].seqid, bed[i].start
    j = k = i
    while bed[j].seqid==seqid and i-j <= 20:
        j -= 1
    j += 1
    while bed[k].seqid==seqid and k-i <= 20:
        k += 1
    k -= 1
    #print j, i, k

    left_dist = abs(bed[j].start - start)
    right_dist = abs(bed[k].start - start)

    query_dr = (max(left_dist, right_dist) / 10000 + 1) * 10000

    # build the csv
    group.insert(0, SyntenyLine("\t".join(str(x) for x in (query, query, "S", k-j, query_dr, "+"))))
    # pad the array to align the columns for each species
    template = ["."] * len(species)

    for g in group:
        prefix = coge_pos[g.anchor][0]
        newname = g.anchor if g.gray=="S" else "-"
        template[species.index(prefix)] = "%s||%d%s" % (newname, g.score, g.gray) 

    gene_names = "\t".join(template)

    # build the link
    url = "http://genomevolution.org/CoGe/GEvo.pl"
    accns = []
    for i, s in enumerate(group):
        j = i+1
        if s.gray=="S":
            accn = "accn%d=%s" % (j, s.anchor)
        else:
            _, dsid, chr, x = coge_pos[s.anchor]
            accn = "dsid%d=%s&chr%d=%s&x%d=%s" % (j, dsid, j, chr, j, x)
        accns.append(accn)
    drups = ["dr%dup=%d" % (i+1, s.dr) for i, s in enumerate(group)]
    drdowns = ["dr%ddown=%d" % (i+1, s.dr) for i, s in enumerate(group)]
    revs = ["rev%d=%d" % (i+1, 1) for i, s in enumerate(group) if s.rev=="-"]
    refs = ["ref%d=%d" % (i+1, 0) for i, s in enumerate(group) if s.anchor!=query]
    extras = ["num_seqs=%d" % len(group), "autogo=1"]
    lyrata_dsid = ["&dsid%d=39219" % (i + 1) for i, s in enumerate(group) \
            if coge_pos[s.anchor][0]=="lyrata" and s.gray=="G"]
    params = "&".join(accns + drups + drdowns + revs + refs + extras + lyrata_dsid)

    return gene_names, url + "?" + params


if __name__ == '__main__':

    # query=>anchors
    rec = collections.defaultdict(list)
    for s in species[1:]:
        attach_species(rec, s)

    fw = file("master_list.tab", "w")
    # header
    print >>fw, "%s\tcoge_link" % ("\t".join(species))
    bed = Bed("bed/all.bed")
    coge_pos = {}
    for x in bed:
        prefix = x.seqid.split(":")[0]
        dsid = DSIDS[prefix]
        seqid = x.seqid.replace(prefix+":", "")
        if prefix in ("athaliana", "grape"): tag = ""
        elif prefix in ("papaya", ): tag = "supercontig_"
        else: tag = "scaffold_"
        coge_pos[x.accn] = (prefix, dsid, tag + seqid, x.start)

    for i, x in enumerate(bed):
        if not x.seqid.startswith("at"): continue
        query = x.accn
        group = list(itertools.chain(*rec[query]))
        csv, gevo_link = build_gevo_link(query, group, i, bed, coge_pos)
        print >>fw, csv + "\t" + gevo_link

