#!/usr/local/bin/python2.6
"""
grab the 'F' rows from the synteny score output,
take the flankers, pad with PAD basepairs, and output
an anchor list for blast_nearby.py

%prog [options] file.synscore

"""
import sys
from bed_utils import Bed
from pyfasta import Fasta

def main(synscore, qbed, qfasta, sbed, sfasta, pad, prefix):
    prefix = prefix.rstrip(".")
    anchor_fh = open("%s.anchors" % prefix, 'w')
    print >>sys.stderr, "writing anchors to %s" % anchor_fh.name
    missing_fh = open("%s.missing" % prefix, 'w')
    print >>sys.stderr, "writing missing sequence to %s" % missing_fh.name
    sfasta = Fasta(sfasta)
    for line in get_anchors(synscore, qbed, sbed, opts.pad):
        # strip the accn and print.
        print >>anchor_fh, "\t".join(map(str, line[1:]))
        fseqid, flank_start, flank_end = line[4:]
        qaccn, f1, f2 = line[0]
        #(f1, f2) == sorted((f1, f2))
        assert flank_start < flank_end, (line)
        seq = str(sfasta[fseqid][flank_start:flank_end]).upper()
        """
        if qaccn == "AT1G01183":
            print >>sys.stderr, "X" * 50
            print >>sys.stderr, qaccn, f1, f2, fseqid, flank_start, flank_end, seq.count('N')
            print >>sys.stderr, "\t".join(map(str, line[1:]))
        """
        print >>missing_fh, "\t".join(map(str, (qaccn, f1, f2, seq.count('N'), seq.count('X'))))


def get_anchors(synscore, qbed, sbed, pad):
    qorder = qbed.get_order()
    sorder = sbed.get_order()

    for line in open(synscore):
        if line[0] == "#": continue
        line = line.rstrip().split("\t")
        if line[2] != 'F': continue
        qaccn, flanker1, flanker2 = line[0], line[1], line[6]

        # returns index, feature. just want feature.
        try:
            q = qorder[qaccn][1]
            f1 = sorder[flanker1][1]
            f2 = sorder[flanker2][1]
        except KeyError, e:
            #print >>sys.stderr, str(e)
            # query and subject gray genes are in the same file.
            # this must be a subject gray gene which we skip.
            continue
        if f1.start > f2.start:
            f1, f2 = f2, f1
        flank_start = max(1, f1.start - pad)
        flank_end = f2.end + pad
        assert f1.seqid == f2.seqid

        yield (qaccn, flanker1, flanker2), q.seqid, q.start, q.end, f1.seqid, flank_start, flank_end


if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    p.add_option('--pad', dest='pad', type='int', default=15000,
                 help="pad region by this many basepairs on either end")
    p.add_option('--qbed', dest='qbed',
                 help="path to query bed file")
    p.add_option('--qfasta', dest='qfasta',
                 help="path to query fasta file")
    p.add_option('--sbed', dest='sbed',
                 help="path to subject bed file")
    p.add_option('--sfasta', dest='sfasta',
                 help="path to subject fasta file")
    p.add_option("--out-prefix", dest='out', help="prefix of"
                 "where to send .anchors and .missing files")
    opts, args = p.parse_args()

    if not (len(args) == 1 and all((opts.qbed, opts.sbed, opts.qfasta, opts.sfasta))):
        sys.exit(not p.print_help())

    qbed = Bed(opts.qbed)
    sbed = Bed(opts.sbed)
    main(args[0], qbed, opts.qfasta, sbed, opts.sfasta, opts.pad, opts.out)
