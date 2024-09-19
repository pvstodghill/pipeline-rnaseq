PIPELINE := $(realpath $(dir $(firstword $(MAKEFILE_LIST))))

export PATH := $(PIPELINE)/howto:${PATH}
export PATH := $(PIPELINE)/stubs:${PATH}

.PHONY:: default all

default: all

include config.mk

########################################################################
# Grab input files
########################################################################

# all:: ${DATA}/inputs/genome.fna
${DATA}/inputs/genome.fna : ${GENOME_FNA}
	@mkdir -p $$(dirname $@)
	cp -a $< $@

# all:: ${DATA}/inputs/genome.gtf
${DATA}/inputs/genome.gtf : ${GENOME_GTF}
	@mkdir -p $$(dirname $@)
	cp -a $< $@

# all:: ${DATA}/inputs/additional.gtf
ifneq (${ADDITIONAL_GTF},)
${DATA}/inputs/additional.gtf : ${ADDITIONAL_GTF}
	@mkdir -p $$(dirname $@)
	cp -a $< $@
else
${DATA}/inputs/additional.gtf :
	@mkdir -p $$(dirname $@)
	touch $@
endif


# all:: $(foreach NAME, $(SAMPLE_NAMES),\
#	   $(foreach Rx, R1 R2,\
#		${DATA}/inputs/${NAME}_${Rx}.fastq.gz))

define rule_copy_fq
${DATA}/inputs/$(1)_$(2).fastq.gz : ${$(2)_$(1)}
	@mkdir -p $$$$(dirname $$@)
	cp -a ${$(2)_$(1)} ${DATA}/inputs/$(1)_$(2).fastq.gz 
endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(foreach Rx, R1 R2,\
	$(eval $(call rule_copy_fq,${NAME},${Rx}))))


########################################################################
# Falco
########################################################################

.PHONY:: falco
all:: falco
falco: $(foreach NAME, $(SAMPLE_NAMES),\
	   $(foreach Rx, R1 R2,\
		${DATA}/falco/${NAME}_${Rx}/summary.txt))

define rule_run_falco
${DATA}/falco/$(1)_$(2)/summary.txt : ${DATA}/inputs/$(1)_$(2).fastq.gz
	D=${DATA}/falco/$$$$(basename $$< .fastq.gz) ; \
	rm -rf $$$$D ; mkdir -p $$$$D ; \
	falco -q -o $$$$D $$<
endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(foreach Rx, R1 R2,\
	$(eval $(call rule_run_falco,${NAME},${Rx}))))

########################################################################
# Fastp
########################################################################

# all:: $(foreach NAME, $(SAMPLE_NAMES),\
#	   $(foreach Rx, R1 R2,\
#		${DATA}/fastp/${NAME}_${Rx}.fastq.gz))

FASTP_THREADS = $(shell if [ ${THREADS} -gt 16 ] ; then echo 16 ; else echo ${THREADS} ; fi)

define rule_run_fastp
${DATA}/fastp/$(1)_R1.fastq.gz ${DATA}/fastp/$(1)_R2.fastq.gz : \
		${DATA}/inputs/$(1)_R1.fastq.gz \
		${DATA}/inputs/$(1)_R2.fastq.gz
	@mkdir -p ${DATA}/fastp
	fastp \
	    --thread ${FASTP_THREADS} \
	     ${FASTP_ARGS} \
	     --in1 ${DATA}/inputs/$(1)_R1.fastq.gz \
	     --in2 ${DATA}/inputs/$(1)_R2.fastq.gz \
	     --json ${DATA}/fastp/$(1).json \
	     --html ${DATA}/fastp/$(1).html \
	     --out1 ${DATA}/fastp/$(1)_R1.fastq.gz \
	     --out2 ${DATA}/fastp/$(1)_R2.fastq.gz \
	     --unpaired1 ${DATA}/fastp/$(1)_u.fastq.gz \
	     --unpaired2 ${DATA}/fastp/$(1)_u.fastq.gz
endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(eval $(call rule_run_fastp,${NAME})))

########################################################################
# Bowtie2
########################################################################

# all:: ${DATA}/bowtie2/genome+phix.fna
${DATA}/bowtie2/genome+phix.fna : ${DATA}/inputs/genome.fna ${PIPELINE}/inputs/phix.fna
	@mkdir -p $(dir $@)
	cat $^ > $@
	bowtie2-build -q --threads ${THREADS} $@ $@

# all:: $(foreach NAME, $(SAMPLE_NAMES), ${DATA}/bowtie2/${NAME}.bam)

define rule_run_bowtie2
${DATA}/bowtie2/$(1).bam : \
		${DATA}/fastp/$(1)_R1.fastq.gz ${DATA}/fastp/$(1)_R2.fastq.gz \
		${DATA}/bowtie2/genome+phix.fna
	@mkdir -p $$$$(dirname $$@)
	bowtie2 --seed 1 --threads ${THREADS} --end-to-end \
		-x ${DATA}/bowtie2/genome+phix.fna \
		-1 ${DATA}/fastp/$(1)_R1.fastq.gz \
		-2 ${DATA}/fastp/$(1)_R2.fastq.gz \
		-S ${DATA}/bowtie2/$(1).sam
	samtools sort -@ ${THREADS} ${DATA}/bowtie2/$(1).sam -o ${DATA}/bowtie2/$(1).bam
	samtools index -@ ${THREADS} ${DATA}/bowtie2/$(1).bam
	rm -f ${DATA}/bowtie2/$(1).sam

endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(eval $(call rule_run_bowtie2,${NAME})))

########################################################################
# Strandedness test
########################################################################

# all:: ${DATA}/strand/targets.gtf
${DATA}/strand/targets.gtf : ${DATA}/inputs/genome.gtf # ${DATA}/inputs/additional.gtf
	@mkdir -p $$(dirname $@)
	cat $^ \
	    | ( fgrep "	gene	" ; true ) \
	    | perl ${PIPELINE}/scripts/sanitize-gtf-for-featureCounts \
	       > $@

A_SAMPLE_NAME = $(firstword ${SAMPLE_NAMES})

${DATA}/strand/forward.txt : ${DATA}/strand/targets.gtf ${DATA}/bowtie2/${A_SAMPLE_NAME}.bam
	@mkdir -p $$(dirname $@)
	featureCounts \
	    -O -p --countReadPairs -B -P -C \
	    -t gene -g gene_id -f \
	    -s 1 \
	    -a ${DATA}/strand/targets.gtf \
	    -o $@ \
	    ${DATA}/bowtie2/${A_SAMPLE_NAME}.bam

${DATA}/strand/reverse.txt : ${DATA}/strand/targets.gtf ${DATA}/bowtie2/${A_SAMPLE_NAME}.bam
	@mkdir -p $$(dirname $@)
	featureCounts \
	    -O -p --countReadPairs -B -P -C \
	    -t gene -g gene_id -f \
	    -s 2 \
	    -a ${DATA}/strand/targets.gtf \
	    -o $@ \
	    ${DATA}/bowtie2/${A_SAMPLE_NAME}.bam

# all:: ${DATA}/strand/results.sh

ifneq (${ORIENTATION},)

${DATA}/strand/results.sh :
	@mkdir -p $$(dirname $@)
	echo ORIENTATION=${ORIENTATION} > $@

else

${DATA}/strand/results.sh : ${DATA}/strand/forward.txt ${DATA}/strand/reverse.txt
	@mkdir -p $$(dirname $@)
	${PIPELINE}/scripts/rnaseq-strand-analysis.pl -o $@ $^

endif

########################################################################
# Make profiles
########################################################################

.PHONY:: profiles
all:: profiles

profiles:: $(foreach NAME, $(SAMPLE_NAMES), ${DATA}/profiles/.done.${NAME})

define rule_make_profiles
${DATA}/profiles/.done.$(1) : \
		${DATA}/bowtie2/$(1).bam \
		${DATA}/strand/results.sh
	@mkdir -p ${DATA}/profiles
	( \
	    . ${DATA}/strand/results.sh ; \
	    case X"$$$$ORIENTATION"X in \
		XforwardX) opt_r= ;; \
		XreverseX) opt_r=-r ;; \
		X*X) echo 1>&2 cannot happen ; exit 1 ; \
	    esac ; \
	    samtools view -h ${DATA}/bowtie2/$(1).bam \
		| ${PIPELINE}/scripts/sam2profiles -2 $$$$opt_r \
						     -e -s -n \
						     -d ${DATA}/profiles \
						     -t $(1) \
	)
	touch ${DATA}/profiles/.done.$(1)
endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(eval $(call rule_make_profiles,${NAME})))

########################################################################
# Make count table
########################################################################

# all:: ${DATA}/strand/targets.gtf
${DATA}/counts/annotations.gtf : ${DATA}/inputs/genome.gtf ${DATA}/inputs/additional.gtf
	@mkdir -p $$(dirname $@)
	cat $^ \
	    | ( fgrep "	${COUNT_FEATURE}	" ; true ) \
	    | perl ${PIPELINE}/scripts/sanitize-gtf-for-featureCounts \
	       > $@

.PHONY:: counts
all:: counts

counts:: $(foreach NAME, $(SAMPLE_NAMES), ${DATA}/counts/counts_${NAME}.txt)

define rule_make_counts
${DATA}/counts/counts_$(1).txt : ${DATA}/bowtie2/$(1).bam \
		${DATA}/counts/annotations.gtf ${DATA}/strand/results.sh
	@mkdir -p ${DATA}/counts
	( \
	    . ${DATA}/strand/results.sh ; \
	    case X"$$$$ORIENTATION"X in \
		XforwardX) s=1 ;; \
		XreverseX) s=2 ;; \
		X*X) echo 1>&2 cannot happen ; exit 1 ; \
	    esac ; \
	    featureCounts \
		-O -p --countReadPairs -B -P -C -s $$$$s \
		-t ${COUNT_FEATURE} -g gene_id -f \
		-s 1 \
		-a ${DATA}/counts/annotations.gtf \
		-o ${DATA}/counts/counts_$(1).txt \
		${DATA}/bowtie2/$(1).bam \
	)
endef

$(foreach NAME, $(SAMPLE_NAMES),\
    $(eval $(call rule_make_counts,${NAME})))

#########################################################################
# Run DESeq2
#########################################################################

${DATA}/deseq2/temp/regions.gtf: ${DATA}/counts/annotations.gtf
	@mkdir -p $$(dirname $@)
	cp $^ $@

all:: $(foreach TAG, ${EXP_NAMES}, ${DATA}/deseq2/results_${TAG}.txt)
all:: $(foreach TAG, ${EXP_NAMES}, ${DATA}/deseq2/results_${TAG}.gff)
all:: $(foreach TAG, ${EXP_NAMES}, ${DATA}/deseq2/changed_${TAG}.gff)


define rule_run_deseq2

${DATA}/deseq2/temp/counts_$(1).txt : ${DATA}/deseq2/temp/regions.gtf \
		$$(foreach NAME, $${CONTROL_SAMPLES_$(1)}, ${DATA}/counts/counts_$${NAME}.txt) \
		$$(foreach NAME, $${TREATMENT_SAMPLES_$(1)}, ${DATA}/counts/counts_$${NAME}.txt)
	@mkdir -p $$$$(dirname $$@)
	$${PIPELINE}/scripts/make-counts-table-from-featurecounts \
	    ${DATA}/deseq2/temp/regions.gtf \
		$$(foreach NAME, $${CONTROL_SAMPLES_$(1)}, $${NAME}:${DATA}/counts/counts_$${NAME}.txt) \
		$$(foreach NAME, $${TREATMENT_SAMPLES_$(1)}, $${NAME}:${DATA}/counts/counts_$${NAME}.txt) \
		> ${DATA}/deseq2/temp/counts_$(1).txt


${DATA}/deseq2/temp/params_$(1).R  : ${DATA}/deseq2/temp/counts_$(1).txt
	@mkdir -p $$$$(dirname $$@)
	$${PIPELINE}/scripts/prep-deseq2 -x -s $${PIPELINE}/scripts \
	    -F parametric \
	    $$(if $${FDR_CUTOFF_$(1)},-f $${FDR_CUTOFF_$(1)}) \
	    $$(if $${PVALUE_CUTOFF_$(1)},-p $${PVALUE_CUTOFF_$(1)}) \
	    -d $$$$(dirname $$@) \
	    -t _$(1) \
	    -c ${DATA}/deseq2/temp/counts_$(1).txt \
	    $$(foreach NAME, $${CONTROL_SAMPLES_$(1)}, $${CONTROL_NAME_$(1)}:$${NAME}) \
	    $$(foreach NAME, $${TREATMENT_SAMPLES_$(1)}, $${TREATMENT_NAME_$(1)}:$${NAME})

${DATA}/deseq2/temp/output_$(1).txt \
${DATA}/deseq2/temp/output-extended_$(1).txt \
: ${DATA}/deseq2/temp/params_$(1).R
	Rscript $${PIPELINE}/scripts/run-deseq2 ${DATA}/deseq2/temp/params_$(1).R


${DATA}/deseq2/results_$(1).txt : ${DATA}/deseq2/temp/regions.gtf \
		${DATA}/deseq2/temp/output-extended_$(1).txt
	cat ${DATA}/deseq2/temp/output-extended_$(1).txt \
	    | $${PIPELINE}/scripts/deseq-output2 -t \
		${DATA}/deseq2/temp/regions.gtf \
		${ALIASES_TXT} \
		> ${DATA}/deseq2/results_$(1).txt


${DATA}/deseq2/results_$(1).gff : ${DATA}/deseq2/temp/output_$(1).txt \
		${DATA}/deseq2/temp/regions.gtf
	cat ${DATA}/deseq2/temp/output_$(1).txt \
	    | $${PIPELINE}/scripts/deseq-output2 -g \
		$$(if $${FOLDCHANGE_CUTOFF_$(1)}, -c $${FOLDCHANGE_CUTOFF_$(1)}) \
		${DATA}/deseq2/temp/regions.gtf \
		${ALIASES_TXT}  \
		> ${DATA}/deseq2/results_$(1).gff

${DATA}/deseq2/changed_$(1).gff : ${DATA}/deseq2/results_$(1).gff
	cat ${DATA}/deseq2/results_$(1).gff \
	    | ( egrep '; colour [23];' ; true) \
	    > ${DATA}/deseq2/changed_$(1).gff

endef

$(foreach TAG, ${EXP_NAMES},\
    $(eval $(call rule_run_deseq2,${TAG})))

########################################################################
# Software versions
########################################################################

all:: ${DATA}/versions.txt

print_version = ( echo + $(1) ; $(1) ; echo ) 2>&1 | tee -a ${DATA}/versions.txt

${DATA}/versions.txt :
	@rm -f ${DATA}/versions.txt
	-@$(call print_version, bowtie2 --version)
	-@$(call print_version, Rscript --version)
	-@$(call print_version, ${PIPELINE}/scripts/print-r-library-versions DESeq2)
	-@$(call print_version, falco --version)
	-@$(call print_version, fastp --version)
	-@$(call print_version, samtools --version | head -n1)
	-@$(call print_version, featureCounts -v)


########################################################################
# Git status
########################################################################

all:: ${DATA}/git-info.txt

${DATA}/git-info.txt :
	-@( \
	    cd ${PIPELINE} ; \
	    echo ; \
	    ( set -x ; git status ) ; \
	    echo ; \
	    ( set -x ; git log -n1 ) \
	) 2>&1 | tee ${DATA}/git-info.txt


########################################################################
# Done!
########################################################################

