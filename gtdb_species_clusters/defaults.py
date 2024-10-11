# criteria for QC'ing genomes
QC_MIN_COMP = 50
QC_MAX_CONT = 10
QC_MIN_QUALITY = 50
QC_SH_EXCEPTION = 80
QC_MIN_PERC_MARKERS = 40
QC_MAX_CONTIGS = 2000
QC_MIN_N50 = 5000
QC_MAX_AMBIGUOUS = 100000

# criteria for defining ANI-based species clusters
ANI_SP = 95.0
AF_SP = 50.0

# criteria used for defining synonyms
ANI_SYNONYMS = 97.0

# minimum ANI estimate from skani sketch estimate for genome
# pair to be fully processsed
SKANI_PREFILTER_THRESHOLD =85.0

# skani presets selecting for different compression factors (k-mer sampling rates)
# None = default which sets c to 125
SKANI_PRESET = None
