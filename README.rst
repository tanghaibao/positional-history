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
