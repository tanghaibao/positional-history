#!/usr/local/bin/python
"""
take the file from synteny_score and the blast run by blast_nearby on the 'F'
genes and tabulate the hits in each region for the gray gene.

%prog [options] file.synscore

"""
import sys
from bed_utils import Bed, BlastLine
import collections
from scipy.spatial.ckdtree import cKDTree
import numpy as np
"""
AT1G01010   GSVIVT01027464001   G   10  440000  +   GSVIVT01027464001
AT1G01010   GSVIVT01019699001   G   4   220000  -   GSVIVT01019699001
AT1G01020   GSVIVT01027464001   S   11  440000  +   GSVIVT01027464001
"""

class SynLine(object):
    __slots__ = ('gray', 'f1', 'tag', 'score', 'width', 'strand', 'f2')
    def __init__(self, line):
        line = line.rstrip().split("\t")
        self.gray = line[0]
        self.f1 = line[1]
        self.tag = line[2]
        self.score = line[3]
        self.width = line[4]
        self.strand = line[5]
        if len(line) > 5:
            self.f2 = line[6]

    def __str__(self):
        return "\t".join(getattr(self, a) for a in self.__slots__[:7])

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.gray)

def get_missing_tag(sl, missing):
    of1, of2 = sorted((sl.f1, sl.f2))
    missing_key = sl.gray, of1, of2
    smiss = missing.get(missing_key)
    tag = ""
    if smiss is None: return tag
    if smiss[0] != 0:
        tag += "N%i" % smiss[0]
    if smiss[1] != 0:
        tag += "X%i" % smiss[1]
    return tag

def main(synscore, qbed, sbed, blasts, trees, tree_dist, missing, blast_tag='B'):
    qorder = qbed.get_order()
    sorder = sbed.get_order()

    for sl in (SynLine(l) for l in open(synscore) if not 'na\t' in l \
               and l.split("\t")[0] in qorder):
        # dont care unless it's an actual flanked gene.
        q, s1, s2 = qorder[sl.gray][1], sorder[sl.f1][1], sorder[sl.f2][1]
        if sl.tag[0] != 'F':
            print str(sl)
            continue
        missing_tag = get_missing_tag(sl, missing)

        try:
            tree = trees[(q.seqid, s1.seqid)]
        except KeyError:
            sl.tag += missing_tag
            print str(sl)
            continue
        smin = min(s1.start, s2.start)
        smax = min(s1.end, s2.end)
        pts = [(q.start, smin), (q.start, smax), (q.end, smax), (q.end, smin)]
        # currently, we only really need to find 1, but in the future, may want
        # to report the blast info for each nearby hit.
        dists, idxs = tree.query(pts, k=2, distance_upper_bound=tree_dist)
        near = []
        for dist, idx in zip(dists, idxs):
            idx = idx[dist != np.inf]
            if len(idx) != 0:
                blast = blasts[(q.seqid, s1.seqid)]
                idx %= float(len(blast))
                #print idx, len(blast)
                idx = np.unique(idx)
                #print idx, len(blast)
                near.extend((blast[i] for i in idx))
                #raw_input("...\n")
                break
        for n in sorted(set(near)):
            sl.tag += blast_tag + str(int(b.score))
        #if near:
            #sl.tag += (blast_tag * len(near))
        sl.tag += missing_tag
        print str(sl)

def read_missing(fmissing):
    d = {}
    for line in open(fmissing):
        line = line.split("\t")
        q, f1, f2 = line[:3]
        f1, f2 = sorted((f1, f2))
        key = (q, f1, f2)
        # q, f1, f2  = N count, X count
        #if q == "AT1G01630":
            #    print >>sys.stderr, key, line[3:]
        d[key] = map(int, line[3:])
    return d

if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser()
    p.add_option('--qbed', dest='qbed',
                 help="path to query bed file")
    p.add_option('--sbed', dest='sbed',
                 help="path to subject bed file")

    p.add_option("--blast", dest="blast",
        help="path to blast file from regions made by extract_nearby.py")
    p.add_option("--dist", dest="dist", type='int',
                 help="distance to check for hits")
    p.add_option("--missing", dest="missing",
                 help="file containing counts of 'N', 'X' for flanking region")

    opts, args = p.parse_args()

    if not all((opts.qbed, opts.sbed, opts.blast, opts.dist)):
        sys.exit(p.print_help())
    if not len(args) == 1:
        sys.exit(p.print_help())

    synscore = args[0]

    missing = read_missing(opts.missing)


    qbed = Bed(opts.qbed)
    sbed = Bed(opts.sbed)

    # blast query, subject are chromosomes.
    blasts = collections.defaultdict(list)
    for line in open(opts.blast):
        b = BlastLine(line)
        key = b.query, b.subject
        blasts[key].append(b)
    trees = {}
    for k, blast_list in blasts.iteritems():
        # since kdtree is points, just add all 4 corners.
        # then in query, can regain the index to blasts by % len(blasts[k])
        xys  = [[b.qstart, b.sstart] for b in blast_list]
        xys += [[b.qstart, b.sstop ] for b in blast_list]
        xys += [[b.qstop,  b.sstart] for b in blast_list]
        xys += [[b.qstop,  b.sstop ] for b in blast_list]
        xys = np.array(xys)
        trees[k] = cKDTree(xys)
    main(synscore, qbed, sbed, blasts, trees, opts.dist, missing)
