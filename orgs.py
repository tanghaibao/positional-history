'''
this will decide what species go into maggies final master spreadsheet
- check DSIDS using CoGe's orgview
- add quota with respect to the expected ploidy level (for example, if poplar
  has its own tetraploidy after athaliana-poplar split, make it 2 - remember the
  quota is relative to athaliana)
'''

from collections import namedtuple

# Represents the configuration of a species, and additional info
# check `dsid` and `tag` using CoGe OrganismView
Species = namedtuple('Species', 'name dsid quota tag')
orgs = [
        Species('athaliana', 42029, 1, ""),
        Species('cacao', 45781, 1, "super_"),
        Species('castor', 42093, 1, ""),
        Species('lyrata', 39129, 1, "scaffold_"),
        Species('papaya', 43486, 1, "supercontig_"),
        Species('peach', 42478, 1, "scaffold_"),
        Species('poplar', 42119, 2, "scaffold_"),
        Species('grape', 43318, 1, ""),
        ]

orgs = dict((x.name, x) for x in orgs)

# you'd only need to modify the species list to be included in the pipeline 
species = "athaliana lyrata papaya poplar grape".split()

DSIDS = dict((x, orgs[x].dsid) for x in species)
quota = dict((x, orgs[x].quota) for x in species)
tags = dict((x, orgs[x].tag) for x in species)

