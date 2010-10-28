#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
the last script to run to get the master spread sheet that contains the
syntelog list, along with annotation info and gevo links
'''

import csv
import itertools
import collections
from bed_utils import Bed
import sys


# get configuration of species list
from orgs import species, DSIDS, quota, tags 


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


def correct_lyrata(filename='At_genes_shouldbe_F_in_lyrata.txt'):
    # read the genes from the file - these are the At genes that should be
    # marked as F in the spreadsheet
    assert species[1]=='lyrata', 'lyrata should be second species on the " +\
                "list, you have %s' % species

    print >>sys.stderr, "read F-gene list for lyrata correction"
    fp = open(filename)
    should_be_F = set(fp.read().split())
    return should_be_F


def create_master_bed(species):
    print >>sys.stderr, "assemble master bed for %s" % species
    fw = open("bed/all.bed", "w")
    for s in species:
        bed = Bed("bed/%s.bed" % s)
        for b in bed:
            seqid = b.seqid
            if s=="athaliana": 
                seqid = seqid.replace("Chr", "")
            b.seqid = "%s:%s" % (s, seqid)
            print >>fw, b


def read_annotations(filename="TAIR9_functional_descriptions"):
    # read in the a. thaliana annotations
    print >>sys.stderr, "read annotation file", filename
    reader = csv.DictReader(open(filename), delimiter="\t")
    annotations = {}
    for row in reader:
        gene = row["Model_name"]
        gene = gene.rsplit(".", 1)[0]
        annotations[gene] = [row["Type"], row["Short_description"]]

    return annotations


def make_header(species):
    mspecies = []
    for s in species:
        mspecies += ["%s-%d" % (s, i+1) for i in range(quota[s])]

    header = 'athaliana\ttype\tdescription\t%s\tcoge_link' % ('\t'.join(mspecies[1:]))

    return header 


def attach_species(rec, s):
    fp = file("results/athaliana_%s.nearby.synteny_score" % s)
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

    left_dist = abs(bed[j].start - start)
    right_dist = abs(bed[k].start - start)

    query_dr = (max(left_dist, right_dist) / 10000 + 1) * 10000

    # build the csv
    group.insert(0, SyntenyLine("\t".join(str(x) for x in \
            (query, query, "S", k-j, query_dr, "+"))))
    # pad the array to align the columns for each species
    species_count = list(itertools.chain.from_iterable([[x] * quota[x] for x in species])) 
    # species_count will look like ['athaliana', 'lyrata', 'poplar', 'poplar', ...],
    # note how poplar gets two counts, according to quota
    template = ["."] * len(species_count) 

    for g in group:
        prefix = coge_pos[g.anchor][0]
        newname = g.anchor if g.gray=="S" else "-"
        colpos = species_count.index(prefix)
        species_count[colpos] = '.' # mark already used
        #template[colpos] = "%s|%d|%s" % (newname, g.score, g.gray) 
        template[colpos] = newname if g.gray=='S' else g.gray
        
    gene_names = template

    # build the link
    url = "http://genomevolution.com/CoGe/GEvo.pl"
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
    lyrata_dsid = ["dsid%d=%d" % (i + 1, DSIDS["lyrata"]) for i, s in enumerate(group) \
            if coge_pos[s.anchor][0]=="lyrata" and s.gray=="S"]
    params = "&".join(accns + drups + drdowns + revs + refs + extras + lyrata_dsid)

    return gene_names, url + "?" + params


if __name__ == '__main__':
    from datetime import date

    # get all the species.bed into the same file for easy loading
    #create_master_bed(species)
    bed = Bed("bed/all.bed")

    # query=>anchors
    rec = collections.defaultdict(list)
    for s in species[1:]:
        attach_species(rec, s)

    annotations = read_annotations()
    should_be_F = correct_lyrata() 

    fw = file("maggie_list-%s.tab" % date.today(),  "w")
    print >>sys.stderr, "writing to", fw.name
    # header
    print >>fw, make_header(species)
    coge_pos = {}
    for x in bed:
        prefix = x.seqid.split(":")[0]
        dsid = DSIDS[prefix]
        seqid = x.seqid.replace(prefix+":", "")
        tag = tags[prefix] # one of "", "supercontig_", "scaffold_"

        coge_pos[x.accn] = (prefix, dsid, tag + seqid, x.start)

    for i, x in enumerate(bed):
        if not x.seqid.startswith("at"): continue
        query = x.accn
        group = list(itertools.chain(*rec[query]))
        gene_names, gevo_link = build_gevo_link(query, group, i, bed, coge_pos)
        query, rest = gene_names[0], gene_names[1:]
        query = query.split('|')[0]
        # apply the correction on lyrata columns so that it matches maggie's
        # plos genetics findings
        if query in should_be_F: 
            rest[0] = 'F'
        csv = [query] + annotations[query] + rest
        print >>fw, "%s\t%s" % ('\t'.join(csv), gevo_link)

