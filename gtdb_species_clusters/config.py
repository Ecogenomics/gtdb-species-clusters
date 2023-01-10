"""
Configuration information that may need updating between releases. Ideally,
this information would be read from a global configuration file used by
all aspects of the GTDB pipeline, but this does not currently exist.
"""

# directory containing BLAST results between LTP database
# and identified 16S rRNA sequences
LTP_DIR = 'rna_ltp_01_2022'

# file containing LTP taxonomy assignments for 16S rRNA sequences
LTP_RESULTS_FILE = 'ssu.taxonomy.tsv'

# ANI criterion to consider two genomes as representing the same strain
ANI_SAME_STRAIN = 99.0
