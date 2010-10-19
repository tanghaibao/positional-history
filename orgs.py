'''
this will decide what species go into maggies final master spreadsheet
- check DSIDS using CoGe's orgview
- add quota with respect to the expected ploidy level (for example, if poplar
  has its own tetraploidy after athaliana-poplar split, make it 2 - remember the
  quota is relative to athaliana)
'''

species = "athaliana lyrata papaya peach poplar grape".split()
DSIDS = [42029, 39129, 43486, 42478, 42119, 43318]
DSIDS = dict(zip(species, DSIDS))
quota = dict(athaliana=1, grape=1, papaya=1, poplar=2, medicago=2, lyrata=1, peach=1)

