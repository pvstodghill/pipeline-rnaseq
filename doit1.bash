#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 0. Set up
# ------------------------------------------------------------------------

rm -rf data

echo 1>&2 '# Initializing data/...'
mkdir -p data/tmp

INPUTS=data/00_inputs
mkdir -p ${INPUTS}

echo 1>&2 '# Making copies of reference genome...'

echo 1>&2 '##' ${INPUTS}/genome.fna "<-" "${REFERENCES_GENOME}"
cp "${REFERENCES_GENOME}" ${INPUTS}/genome.fna
echo 1>&2 '##' ${INPUTS}/annotation.gff "<-" "${REFERENCES_ANNOTATION_GFF}"
cp "${REFERENCES_ANNOTATION_GFF}" ${INPUTS}/annotation.gff
echo 1>&2 '##' ${INPUTS}/annotation.gtf "<-" "${REFERENCES_ANNOTATION_GTF}"
cp "${REFERENCES_ANNOTATION_GTF}" ${INPUTS}/annotation.gtf

echo 1>&2 '# Making copies of raw reads...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R1.fastq.gz "<-" "${SAMPLES_R1[$i]}"
    cp "${SAMPLES_R1[$i]}" ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R2.fastq.gz "<-" "${SAMPLES_R2[$i]}"
	cp "${SAMPLES_R2[$i]}" ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

if [ "$PACKAGES_FROM" = howto ] ; then
    echo 1>&2 '# Ensuring entries in packages.yaml are downloaded...'
    (
	set -x
	./howto/howto -f packages.yaml -p '*'
    )
fi

# ------------------------------------------------------------------------
# Step 1. FASTQC, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FASTQC on raw reads'

FASTQC1=data/01_fastqc

rm -rf ${FASTQC1}

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${FASTQC1}/${SAMPLES_NAME[$i]}_R1
    mkdir -p ${FASTQC1}/${SAMPLES_NAME[$i]}_R1
    fastqc -t ${THREADS} \
	   -o ${FASTQC1}/${SAMPLES_NAME[$i]}_R1 \
	   ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTQC1}/${SAMPLES_NAME[$i]}_R2
	mkdir -p ${FASTQC1}/${SAMPLES_NAME[$i]}_R2
	fastqc -t ${THREADS} \
	       -o ${FASTQC1}/${SAMPLES_NAME[$i]}_R2 \
	       ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Step 2. Run FASTP on Illumina reads
# ------------------------------------------------------------------------

echo 1>&2 '# Clean-up Illumina reads...'

FASTP=data/02_fastp

rm -rf ${FASTP}
mkdir -p ${FASTP}

for i in $SAMPLES_INDICES ; do
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R2.fastq.gz
	fastp \
		 --thread ${THREADS} \
		 --trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT \
		 --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
		 --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
		 --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
		 --in2 ${INPUTS}/raw_${i}_R2.fastq.gz \
		 --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 --out2 ${FASTP}/trimmed_${i}_R2.fastq.gz \
		 --unpaired1 ${FASTP}/unpaired_$i.fastq.gz \
		 --unpaired2 ${FASTP}/unpaired_$i.fastq.gz
    else
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R1.fastq.gz
	fastp \
		 --thread ${THREADS} \
		 --trim_front1 1 --adapter_sequence CTGTCTCTTATACACATCT \
		 --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
		 --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
		 --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
		 --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Step 3. FASTQC, round 2
# ------------------------------------------------------------------------

if [ "$RUN_FASTQC_AFTER_FASTP" ] ; then

    echo 1>&2 '# Running FASTQC on trimmed reads'

    FASTQC2=data/03_fastqc

    rm -rf ${FASTQC2}

    for i in $SAMPLES_INDICES ; do
	echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
	mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
	fastqc -t ${THREADS} \
	       -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R1 \
	       ${FASTP}/trimmed_${i}_R1.fastq.gz
	if [ "$PE" ] ; then
	    echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	    mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	    fastqc -t ${THREADS} \
		   -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R2 \
		   ${FASTP}/trimmed_${i}_R2.fastq.gz
	fi
    done

fi

# ------------------------------------------------------------------------
# Step 4. Run Bowtie2
# ------------------------------------------------------------------------

echo 1>&2 '# Indexing genome...'

BOWTIE2=data/04_bowtie2
rm -rf ${BOWTIE2}
mkdir -p ${BOWTIE2}

echo 1>&2 '##' ${BOWTIE2}/genome+phix.fna
cat ${INPUTS}/genome.fna inputs/phix.fna > ${BOWTIE2}/genome+phix.fna
bowtie2-build -q --threads ${THREADS} \
	 ${BOWTIE2}/genome+phix.fna \
	 ${BOWTIE2}/genome+phix.fna

echo 1>&2 '# Aligning the reads to the genomes...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${BOWTIE2}/aligned_$i.bam
    echo 1>&2 '###' bowtie2...
    if [ "$PE" ] ; then
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/genome+phix.fna \
		 -1 ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 -2 ${FASTP}/trimmed_${i}_R2.fastq.gz \
		 -S ${BOWTIE2}/aligned_$i.sam
    else
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/genome+phix.fna \
		 -U ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 -S ${BOWTIE2}/aligned_$i.sam
    fi
    echo 1>&2 '###' samtools sort
    samtools sort -@ ${THREADS} ${BOWTIE2}/aligned_$i.sam -o ${BOWTIE2}/aligned_$i.bam
    rm -f ${BOWTIE2}/aligned_$i.sam
    echo 1>&2 '###' samtools index
    samtools index -@ ${THREADS} ${BOWTIE2}/aligned_$i.bam
done

# ------------------------------------------------------------------------
# Step 5. Generating GFF files + Making profiles
# ------------------------------------------------------------------------

PROFILES=data/05_profiles
BAMGFF=${PROFILES}/tmp

rm -rf ${PROFILES}
mkdir -p ${BAMGFF}

echo 1>&2 '# Converting the .bam files to .gff files'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${BAMGFF}/aligned_$i.gff.gz
    if [ "$PE" ] ; then
	samtools view -h ${BOWTIE2}/aligned_$i.bam \
	    | ./scripts/sam2gff -r -p \
	    | gzip -c > ${BAMGFF}/aligned_$i.gff.gz
    else
	samtools view -h ${BOWTIE2}/aligned_$i.bam \
	    | ./scripts/sam2gff -r -s \
	    | gzip -c > ${BAMGFF}/aligned_$i.gff.gz
    fi
done

echo 1>&2 '# Making profiles'

for i in $SAMPLES_INDICES ; do
    name=${SAMPLES_NAME[i]}
    echo 1>&2 '##' $i':' ${PROFILES}/$name'*'.profile

    gzip -dc ${BAMGFF}/aligned_$i.gff.gz \
	| ./scripts/gff2profiles -e -s -n -d ${PROFILES} \
				 ${BOWTIE2}/genome+phix.fna $name
done

# ------------------------------------------------------------------------
# Step 6. Make count tables
# ------------------------------------------------------------------------

echo 1>&2 '# Making count tables'

COUNTS=data/06_counts
rm -rf ${COUNTS}
mkdir -p ${COUNTS}

FEATURECOUNTS_ARGS=
FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

FEATURECOUNTS_ARGS+=" -O" # Assign reads to all their overlapping features

#FEATURECOUNTS_ARGS+=" -M" # all multi-mapping reads reported alignments will be counted
##FEATURECOUNTS_ARGS+=" --fraction" # Assign fractional counts to features
#FEATURECOUNTS_ARGS+=" --primary" # Count primary alignments only !(0x100)

FEATURECOUNTS_ARGS+=" -s 1" # stranded
#FEATURECOUNTS_ARGS+=" -s 2" # reverse-stranded


FEATURECOUNTS_ARGS+=" -p" # fragments (or pairs) will be counted instead of reads (<v2.0.2)
#FEATURECOUNTS_ARGS+=" -p" # input data contains paired-end reads. (>=v2.0.2)
#FEATURECOUNTS_ARGS+=" --countReadPairs" # Count read pairs (fragments) instead of reads

FEATURECOUNTS_ARGS+=" -B" # Only count read pairs that have both ends aligned.
FEATURECOUNTS_ARGS+=" -P" # Check validity of paired-end distance
FEATURECOUNTS_ARGS+=" -C" # Only count concordant reads

FEATURECOUNTS_ARGS+=" -T ${THREADS}"

FEATURECOUNTS_ARGS+=" "
FEATURECOUNTS_ARGS+=" "


# FIXME: cp local/annotations.gtf ${COUNTS}/annotation.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation.gtf

for i in $SAMPLES_INDICES ; do

    featureCounts $FEATURECOUNTS_ARGS \
    		  -a ${COUNTS}/annotation.gtf \
    		  -o ${COUNTS}/counts_$i.txt \
    		  ${BOWTIE2}/aligned_$i.bam

done


# ------------------------------------------------------------------------
# Step 7. Run DESeq2
# ------------------------------------------------------------------------

echo 1>&2 '# Running DESeq2...'

CHANGE_CUTOFF=2.0
TAG=_WT-5255

DESEQ2=data/07_deseq2

rm -rf ${DESEQ2}
mkdir -p ${DESEQ2}/temp

cp ${COUNTS}/annotation.gtf ${DESEQ2}/temp/regions.gtf

./scripts/make-counts-table-from-featurecounts \
    ${DESEQ2}/temp/regions.gtf \
    52551:${COUNTS}/counts_0.txt \
    52552:${COUNTS}/counts_1.txt \
    52553:${COUNTS}/counts_2.txt \
    WT1:${COUNTS}/counts_3.txt \
    WT2:${COUNTS}/counts_4.txt \
    WT3:${COUNTS}/counts_5.txt \
    > ${DESEQ2}/temp/counts.txt

# ------------------------------------------------------------------------

./scripts/prep-deseq2 -x -s ./scripts \
		      -F parametric \
		      -d ${DESEQ2}/temp \
		      -t ${TAG} \
		      -c ${DESEQ2}/temp/counts.txt \
 		     WT:WT1 WT:WT2 WT:WT3 \
 		     5255:52551 5255:52552 5255:52553


# ------------------------------------------------------------------------

Rscript ./scripts/run-deseq2 \
	${DESEQ2}/temp/params${TAG}.R

# ------------------------------------------------------------------------

cat ${DESEQ2}/temp/output-extended${TAG}.txt \
    | ./scripts/deseq-output2results > ${DESEQ2}/results${TAG}.txt

cat ${DESEQ2}/temp/output${TAG}.txt \
    | ./scripts/deseq-output2gff ${ACCESSION} ${CHANGE_CUTOFF} \
				 > ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff
set +e # if egrep matches nothing
cat ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff \
    | egrep '; colour [23];' > ${DESEQ2}/changed${TAG}_${CHANGE_CUTOFF}.gff
set -e

# ------------------------------------------------------------------------
# Step 8. Compute some high-level stats
# ------------------------------------------------------------------------

echo 1>&2 '# Generating statistics...'

STATS=data/08_stats

rm -rf ${STATS}
mkdir -p ${STATS}/temp

cat ${INPUTS}/annotation.gff inputs/phix.gff \
    | egrep -v '^#' | fgrep -v $'\t'region$'\t' \
			    > ${STATS}/temp/regions.gff

./scripts/make-counts-table-from-gffs \
    -t -f -u \
    ${STATS}/temp/regions.gff \
    52551:${BAMGFF}/aligned_0.gff.gz \
    52552:${BAMGFF}/aligned_1.gff.gz \
    52553:${BAMGFF}/aligned_2.gff.gz \
    WT1:${BAMGFF}/aligned_3.gff.gz \
    WT2:${BAMGFF}/aligned_4.gff.gz \
    WT3:${BAMGFF}/aligned_5.gff.gz \
    > ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Package the results
# ------------------------------------------------------------------------

RESULTS=$(date +%Y-%m-%d)-results
echo 1>&2 '# Creating '$RESULTS'...'

rm -rf ${RESULTS}
mkdir -p ${RESULTS}

cp ${PROFILES}/*.profile ${RESULTS}
cp ${DESEQ2}/results_*.txt  ${RESULTS}

for f in ${DESEQ2}/*.gff ; do
    name=$(basename $f .gff)
    cat $f | ./scripts/split-gff -n -d ${RESULTS} $name
done

cp ${STATS}/stats.txt ${RESULTS}/stats.txt
# ./scripts/tsv2xlsx -o ${RESULTS}/stats.xlsx ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
