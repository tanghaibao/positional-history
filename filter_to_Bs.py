import sys
from urllib import urlencode

master_list_tab = sys.argv[1]

def parse_link(idx, link):
    d = dict(x.split("=") for x in link.split("?")[1].split("&"))
    exclude = [x for x in range(2, 6) if x != idx]
    for ex in exclude:
        for k in ('dr%iup', 'dr%idown', 'accn%i', 'chr%i', 'x%i', 'ref%i', 'rev%i', 'dsid%i'):
            try: del d[k % ex]
            except KeyError: pass
    for k in ('dr%iup', 'dr%idown', 'accn%i', 'chr%i', 'x%i', 'ref%i', 'rev%i', 'dsid%i'):
        try:
            d[k % 2] = d[k % idx]
            if idx != 2:
                del d[k % idx]
        except KeyError: pass

    d['num_seqs'] = 2
    d['dr1up'] = 100
    d['dr1down'] = 100
    d['dr2up'] = 20000
    d['dr2down'] = 20000
    d['prog'] = 'blastz'
    include = urlencode(d)
    return "%s?%s" % (link.split("?")[0], include)


for i, line in enumerate(open(sys.argv[1])):
    if i == 0:
        print line,
        continue
    tline = line.split("\t")

    fixer = 0
    for col in (1, 2, 3, 4):
        if "FB" in tline[col] and not "BN" in tline[col]:
            #print line,
            link = tline[-1]
            #print link
            new_link = parse_link(1 + col + fixer, link)
            #print new_link
            tline[-1] = new_link
            print "\t".join(tline)
            break
        elif tline[col] == ".":
            fixer += 1

