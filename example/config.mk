DATA := data

THREADS := $(shell nproc)

GENOME_FNA := example/NC_004632.fna
GENOME_GTF := example/NC_004632.gtf

ADDITIONAL_GTF:=

SAMPLE_NAMES := WT1 WT2 WT3 52551 52552 52553

R1_WT1 := example/WT1_R1.fastq.gz
R2_WT1 := example/WT1_R2.fastq.gz
R1_WT2 := example/WT2_R1.fastq.gz
R2_WT2 := example/WT2_R2.fastq.gz
R1_WT3 := example/WT3_R1.fastq.gz
R2_WT3 := example/WT3_R2.fastq.gz
R1_52551 := example/52551_R1.fastq.gz
R2_52551 := example/52551_R2.fastq.gz
R1_52552 := example/52552_R1.fastq.gz
R2_52552 := example/52552_R2.fastq.gz
R1_52553 := example/52553_R1.fastq.gz
R2_52553 := example/52553_R2.fastq.gz

FASTP_ARGS := --trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT

# ORIENTATION := forward # reads map to "sense"
# ORIENTATION := reverse # reads map to "antisense"
ORIENTATION := # we'll figure it out

COUNT_FEATURE := gene

EXP_NAMES := normal unfiltered

FDR_CUTOFF_normal := 0.05
PVALUE_CUTOFF_normal := -1
FOLDCHANGE_CUTOFF_normal := 2.0
CONTROL_NAME_normal := WT
CONTROL_SAMPLES_normal := WT1 WT2 WT3
TREATMENT_NAME_normal := Delta5255
TREATMENT_SAMPLES_normal := 52551 52552 52553

FDR_CUTOFF_unfiltered := -1
PVALUE_CUTOFF_unfiltered := -1
FOLDCHANGE_CUTOFF_unfiltered := 1.5
CONTROL_NAME_unfiltered := ${CONTROL_NAME_normal}
CONTROL_SAMPLES_unfiltered := ${CONTROL_SAMPLES_normal}
TREATMENT_NAME_unfiltered := ${TREATMENT_NAME_normal}
TREATMENT_SAMPLES_unfiltered := ${TREATMENT_SAMPLES_normal}

