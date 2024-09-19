########################################################################
# Basic parameters
########################################################################

# # Where to put output files
# DATA := data

# # How many threads to use
# THREADS := $(shell nproc)

# # The genome
# GENOME_FNA := fixme.fna
# GENOME_GTF := fixme.gtf
# # Any additional annotations
# ADDITIONAL_GTF:=

########################################################################
# Illumina parameters
########################################################################

# SAMPLE_NAMES := WT1 WT2 WT3 Mutant1 Mutant2 Mutant3

# R1_WT1 := /path/to/raw/reads/WT1_R1.fastq.gz
# R2_WT1 := /path/to/raw/reads/WT1_R2.fastq.gz
# R1_WT2 := /path/to/raw/reads/WT2_R1.fastq.gz
# R2_WT2 := /path/to/raw/reads/WT2_R2.fastq.gz
# R1_WT3 := /path/to/raw/reads/WT3_R1.fastq.gz
# R2_WT3 := /path/to/raw/reads/WT3_R2.fastq.gz
# R1_Mutant1 := /path/to/raw/reads/Mutant1_R1.fastq.gz
# R2_Mutant1 := /path/to/raw/reads/Mutant1_R2.fastq.gz
# R1_Mutant2 := /path/to/raw/reads/Mutant2_R1.fastq.gz
# R2_Mutant2 := /path/to/raw/reads/Mutant2_R2.fastq.gz
# R1_Mutant3 := /path/to/raw/reads/Mutant3_R1.fastq.gz
# R2_Mutant3 := /path/to/raw/reads/Mutant3_R2.fastq.gz

# # FASTP_ARGS := --trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT

# # How do the reads map to the genome?
# # ORIENTATION := forward # reads map to "sense"
# # ORIENTATION := reverse # reads map to "antisense"
# ORIENTATION := # we'll figure it out heuristically

########################################################################
# Diff exp parameters
########################################################################

# COUNT_FEATURE := gene

# EXP_NAMES := normal unfiltered

# FDR_CUTOFF_normal := 0.05
# PVALUE_CUTOFF_normal := -1
# FOLDCHANGE_CUTOFF_normal := 2.0
# CONTROL_NAME_normal := WT
# CONTROL_SAMPLES_normal := WT1 WT2 WT3
# TREATMENT_NAME_normal := Mutant
# TREATMENT_SAMPLES_normal := Mutant1 Mutant2 Mutant3

# FDR_CUTOFF_unfiltered := -1
# PVALUE_CUTOFF_unfiltered := -1
# FOLDCHANGE_CUTOFF_unfiltered := 1.5
# CONTROL_NAME_unfiltered := ${CONTROL_NAME_normal}
# CONTROL_SAMPLES_unfiltered := ${CONTROL_SAMPLES_normal}
# TREATMENT_NAME_unfiltered := ${TREATMENT_NAME_normal}
# TREATMENT_SAMPLES_unfiltered := ${TREATMENT_SAMPLES_normal}

