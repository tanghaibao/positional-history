Pipeline for positional history of A. thaliana genes
=====================================================

This folder contains complete pipeline to run positional history of A. thaliana
genes to generate spreadsheets, which are then manually proofed with the final
data presented on http://biocon.berkeley.edu/athaliana.

The scripts are mostly internal and we do not make attempts to support the
usage for external research. The entire pipeline is run as command lines in the
order specified in ``run.py``.

Modify ``orgs.py`` to select only subset of species to be included in the final
spreadsheet, which will be something like ``master_list-2010-10-21.tab``. Before
running the pipeline, also make sure the following files are in place. 

* ``bed/species.bed``, BED-format (``chr``, ``start``, ``stop``, ``name``)
* ``fastas/species.fa`` for the genome FASTA file, and ``fastas/species.cds``
  for only the CDS FASTA file
* ``blasts/athaliana_species.blastz``, BLAST tabular format (either through
  BLAST -m8 or use LASTZ wrapper `here
  <http://github.com/tanghaibao/bio-pipeline/tree/master/lastz_wrapper/>`_)


Notes
------
``blast_to_raw.py`` has ``strip_names`` in the BLAST parsing routine. Sometimes you
want to strip names for the query but not subject, other times you want to strip
names for both. Modify around Line 80. 

Known issues
-------------
+ RNA genes included in bed, but not blasted.
+ not always getting tightest flankers
+ not sure of gene-set for thaliana.
