# Intron Marker Pipeline (IMP) #
## Sanderson L, Perumal V, Bett K (released 2009) ##

This pipeline has been developed as automated tool for design of primers which flank predicted intron location in Expressed Sequence Tags (ESTs) for species with limited to no genomic sequence data available. All that is required for input is a FASTA file containing the EST sequence and another FASTA file containing any genomic sequence to be used for intron prediction. This allows optimal flexibility for users, allowing them to control the exact sequences used.


**An article detailing IMP and comparing it to current programs available has been submitted. Check back here in the future for a reference.**


---


Please click the Download Tab to download this program. It is an operating system-independent perl script requiring little installation.

Dependencies include:
  1. RepeatMasker -http://www.repeatmasker.org/
  1. BLAT -http://www.kentinformatics.com/products.html
  1. GeneSeqer -http://deepc2.psi.iastate.edu/bioinformatics2go/gs/download.html
  1. Primer3 command-line utility -http://primer3.sourceforge.net/releases.php

Installation Instructions:
  1. Download the current distribution
  1. Unpack (will create it’s own containing folder) using _tar –zxvf filename_ in the terminal
  1. Install ReapeatMasker, BLAT, GeneSeqer and Primer3 command line utilities and note down the path to and including the executable for each one (ex:/usr/local/share/applications/Primer3/src/primers\_core)
  1. Run perl install.pl in the base directory of IMP
  1. Enter the paths already noted down when prompted