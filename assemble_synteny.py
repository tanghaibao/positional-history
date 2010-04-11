#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import csv
import itertools
import collections

"""
import bitly
api = bitly.Api(login="tanghaibao", apikey="R_51b7c08fbb167cc934ca6ae7fd9e8b2d")
def shorten(a, url):
    return a.shorten(url)
"""

quota = dict(grape=1, papaya=1, poplar=2, medicago=2, lyrata=1, peach=1)
#species = "lyrata papaya poplar medicago grape".split()
species = "lyrata papaya peach grape".split()


def coge_url(accns, revs, drs):
    assert len(accns) == len(revs)
    g = ["accn%d=%s" % (i+1, g) for i, g in enumerate(accns)]
    noref_all = ["ref%d=0" % (i+1) for i in xrange(1, len(accns))]
    revs = ["rev%d=1" % (i+1) for i, r in enumerate(revs) if r=="-"]
    drs = ["drup%d=%s&drdown%d=%s" % (i+1, r, i+1, r) for i, r in enumerate(drs)]
    lyrata_dsid = "&dsid%d=39129" % (species.index("lyrata") + 1)
    extra = "&".join(["&".join(noref_all), "&".join(revs), "&".join(drs), lyrata_dsid])
    gevo_url = "http://genomevolution.org/CoGe/GEvo.pl?%s&num_seqs=%d&autogo=1&%s" % \
        ("&".join(g), len(accns), extra)
    return gevo_url

def attach_species(rec, s):
    fp = file("results/athaliana_%s.synteny_score" % s)
    q = quota[s] # get top n hits (ranked by synteny score)
    qrec = collections.defaultdict(list)
    for row in fp:
        query, anchor, dr, gray, orientation, score = row.split()
        if not query.startswith("AT"): continue
        if anchor=="na": anchor, dr, gray, orientation, score = "", "", "", "", ""
        if len(qrec[query]) < q:
            qrec[query].append((anchor, dr, gray, orientation, score))

    for q, v in qrec.iteritems():
        rec[q].append(sorted(v))


if __name__ == '__main__':

    # query=>anchors
    rec = collections.defaultdict(list)
    for i, s in enumerate(species):
        attach_species(rec, s)

    fw = file("master_list.csv", "w")
    # header
    print >>fw, "athaliana,%s,coge_link" % (",".join(species))
    for q, v in sorted(rec.items()):
        if not q.startswith("AT"): continue
        
        all_genes = [[x[0] for x in cv] for cv in v]
        all_genes = [q] + list(itertools.chain(*all_genes))
        
        revs = [[x[3] for x in cv] for cv in v]
        revs = ["+"] + list(itertools.chain(*revs))

        drs =  [[x[1] for x in cv] for cv in v]
        drs = ["75000"] + list(itertools.chain(*drs))

        v = ["|".join("%s(%s%s)" % (a, s, g) for (a, d, g, o, s) in cv) for cv in v]
        #print all_genes
        print >>fw, "%s,%s,%s" % (q, ",".join(v), coge_url(all_genes, revs, drs)) 

    fw.close()
