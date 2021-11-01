#! /bin/bash

. doit-preamble.bash

BOWTIE2=data/04_bowtie2

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
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
