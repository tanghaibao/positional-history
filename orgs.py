'''
this will decide what species go into maggies final master spreadsheet
- check DSIDS using CoGe's orgview
- add quota with respect to the expected ploidy level (for example, if poplar
  has its own tetraploidy after athaliana-poplar split, make it 2 - remember the
  quota is relative to athaliana)
'''

from collections import namedtuple

# Represents the configuration of a species, and additional info
Species = namedtuple('Species', 'name dsid quota')
orgs = [
        Species('athaliana', 42029, 1),
        Species('cacao', 45781, 1),
        Species('lyrata', 39129, 1),
        Species('papaya', 43486, 1),
        Species('peach', 42478, 1),
        Species('poplar', 42119, 2),
        Species('grape', 43318, 1),
        ]

orgs = dict((x.name, x) for x in orgs)

# you'd only need to modify the species list to be included in the pipeline 
species = "athaliana lyrata papaya cacao peach grape".split()

DSIDS = dict((x, orgs[x].dsid) for x in species)
quota = dict((x, orgs[x].quota) for x in species)

