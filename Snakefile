import os
import glob

configfile: "config.yaml"

SAMPLES_NAMES = config['samples'].keys()
EXPS_NAMES = config['exps'].keys()

DATA=config['data'] if 'data' in config else "data"
PIPELINE=os.path.dirname(workflow.snakefile)

# ------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------

FALCO_RESULTS = \
    expand(DATA+"/falco/{name}_{rx}",
           name = SAMPLES_NAMES, rx=['R1','R2']),

PROFILE_RESULTS = \
    expand(DATA+"/profiles/.done.{name}",
           name = SAMPLES_NAMES),

SOUS_RESULTS = DATA+"/sous/.done"

DESEQ2_RESULTS = \
    expand(DATA+"/deseq2/changed_{tag}.gff", tag=EXPS_NAMES) \
    + expand(DATA+"/deseq2/results_{tag}.gff", tag=EXPS_NAMES) \
    + expand(DATA+"/deseq2/results_{tag}.txt", tag=EXPS_NAMES)

VERSIONS_RESULTS = \
    expand(DATA+"/versions/{pkg}.txt",
           pkg = list(map((lambda p: os.path.splitext(os.path.basename(p))[0] ),glob.glob(PIPELINE+"/envs/*.yaml")))),

GIT_RESULTS = DATA+"/git.log"

rule all:
    input:
        FALCO_RESULTS, \
        DESEQ2_RESULTS, \
        PROFILE_RESULTS, \
        SOUS_RESULTS, \
        VERSIONS_RESULTS, \
        GIT_RESULTS,

        
# ------------------------------------------------------------------------
# Collect inputs
# ------------------------------------------------------------------------

rule make_inputs:
    input:
        DATA+"/inputs/genome.fna",
        DATA+"/inputs/genome.gtf",
        DATA+"/inputs/additional.gtf" if 'additional_gtf' in config['genome'] else [],
        expand(DATA+"/inputs/{name}_{rx}.fastq.gz", name = SAMPLES_NAMES, rx = ['R1','R2'])

if 'gbk' in config['genome']:
    rule copy_genome_gbk:
        input: config['genome']['gbk']
        output: DATA+"/inputs/genome.gbk"
        shell: "cat {input} > {output}"
    rule make_genome_fna:
        input: DATA+"/inputs/genome.gbk"
        output: DATA+"/inputs/genome.fna"
        shell:
            """
            {PIPELINE}/scripts/gbk2fna < {input} > {output}
            """
    rule make_genome_gff:
        input: DATA+"/inputs/genome.gbk"
        output: DATA+"/inputs/genome.gff"
        shell:
            """
            {PIPELINE}/scripts/gbk2gff < {input} > {output}
            """
    rule make_genome_gtf:
        input: DATA+"/inputs/genome.gff"
        output: DATA+"/inputs/genome.gtf"
        shell:
            """
            cat {input} \
            | fgrep $'\t'gene$'\t' \
            | {PIPELINE}/scripts/gff2gtf \
            > {output}
            """

if 'fna' in config['genome']:
    rule copy_genome_fna:
        input: config['genome']['fna']
        output: DATA+"/inputs/genome.fna"
        shell: "cat {input} > {output}"

if 'gtf' in config['genome']:
    rule copy_genome_gtf:
        input: config['genome']['gtf']
        output: DATA+"/inputs/genome.gtf"
        shell: "cat {input} > {output}"

if 'additional_gtf' in config['genome']:
    rule copy_additional_gtf:
        input: config['genome']['additional_gtf']
        output: DATA+"/inputs/additional.gtf"
        shell: "cat {input} > {output}"
        
rule make_annotations_gtf:
    input: 
        DATA+"/inputs/genome.gtf",
        DATA+"/inputs/additional.gtf" if 'additional_gtf' in config['genome'] else [],
    output: DATA+"/inputs/annotations.gtf"
    shell: "cat {input} > {output}"

SAMPLES_PREFIX=os.path.expanduser(config['samples_dir']+'/') if 'samples_dir' in config else ''

rule copy_fq:
    input: lambda wildcards: SAMPLES_PREFIX+config['samples'][wildcards.name][wildcards.rx]
    output: DATA+"/inputs/{name}_{rx}.fastq.gz"
    shell: "cat {input} > {output}"
    
# ------------------------------------------------------------------------
# Run Falco
# ------------------------------------------------------------------------

rule run_falco:
    input: FALCO_RESULTS

rule falco_version:
    output: DATA+"/versions/falco.txt"
    conda: "envs/falco.yaml"
    shell:
        """
        falco --version 2>&1 | tee {output}
        """    

rule run_falco_once:
    input: DATA+"/inputs/{name}_{rx}.fastq.gz"
    output: directory(DATA+"/falco/{name}_{rx}")
    threads: 1
    conda: "envs/falco.yaml"
    shell:
        """
        falco -q -o {output} {input}
        """

# ------------------------------------------------------------------------
# Run Fastp
# ------------------------------------------------------------------------

rule run_fastp:
    input: expand(DATA+"/fastp/{name}_{rx}.fastq.gz", \
                  name = SAMPLES_NAMES, \
                  rx = ['R1','R2'])

rule fastp_version:
    output: DATA+"/versions/fastp.txt"
    conda: "envs/fastp.yaml"
    shell:
        """
        fastp --version 2>&1 | tee {output}
        """    

rule run_fastp_one:
    input:
        R1=DATA+"/inputs/{name}_R1.fastq.gz",
        R2=DATA+"/inputs/{name}_R2.fastq.gz"
    output:
        R1=DATA+"/fastp/{name}_R1.fastq.gz",
        R2=DATA+"/fastp/{name}_R2.fastq.gz",
        unpaired=DATA+"/fastp/{name}_u.fastq.gz",
        json=DATA+"/fastp/{name}.json",
        html=DATA+"/fastp/{name}.html"
    params: config['fastp_adapter_args'] if 'fastp_adapter_args' in config else ""
    threads: 16
    conda: "envs/fastp.yaml"
    shell:
        """
        fastp \
            --thread {threads} \
            {params} \
            --in1 {input.R1} --in2 {input.R2} \
            --json {output.json} --html {output.html} \
            --out1 {output.R1} --out2 {output.R2} \
            --unpaired1 {output.unpaired} --unpaired2 {output.unpaired}
        """

# ------------------------------------------------------------------------
# Run Bowtie2
# ------------------------------------------------------------------------

rule run_bowtie2:
    input: expand(DATA+"/bowtie2/{name}.sam", \
                  name = SAMPLES_NAMES)

rule bowtie2_version:
    output: DATA+"/versions/bowtie2.txt"
    conda: "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 --version 2>&1 | tee {output}
        """    

rule make_genome_phix:
    input:
        DATA+"/inputs/genome.fna",
        PIPELINE+"/inputs/phix.fna"
    output:
        DATA+"/bowtie2/genome+phix.fna"
    threads: 9999
    conda: "envs/bowtie2.yaml"
    shell:
        """
        cat {input} > {output}
        bowtie2-build -q --threads {threads} {output} {output}
        """

rule run_bowtie2_once:
    input:
        fna=DATA+"/bowtie2/genome+phix.fna",
        R1=DATA+"/fastp/{name}_R1.fastq.gz",
        R2=DATA+"/fastp/{name}_R2.fastq.gz",
    output:
        sam=temp(DATA+"/bowtie2/{name}.sam"),
    threads: 9999
    conda: "envs/bowtie2.yaml"
    shell:
        """
        bowtie2 --seed 1 --threads {threads} --end-to-end \
                -x {input.fna} -1 {input.R1} -2 {input.R2} \
                 -S {output.sam}
        """

# ------------------------------------------------------------------------
# Run Samtools: .sam -> .bam
# ------------------------------------------------------------------------

rule run_samtools:
    input: expand(DATA+"/bowtie2/{name}.bam", \
                  name = SAMPLES_NAMES)

rule samtools_version:
    output: DATA+"/versions/samtools.txt"
    conda: "envs/samtools.yaml"
    shell:
        """
        samtools --version 2>&1 | tee {output}
        """    

rule make_bam:
    input:
        sam=DATA+"/bowtie2/{name}.sam"
    output:
        bam=DATA+"/bowtie2/{name}.bam"
    threads: 9999
    conda: "envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        samtools index -@ {threads} {output}
        """

# ------------------------------------------------------------------------
# Test strandedness
# ------------------------------------------------------------------------

rule test_strand:
    input: DATA+"/strand/results.sh"

rule make_strand_targets:
    input: DATA+"/inputs/annotations.gtf"
    output: DATA+"/strand/targets.gtf"
    shell:
        """
        cat {input} \
            | perl {PIPELINE}/scripts/sanitize-gtf-for-featureCounts \
        	-f gene \
                > {output}
        """

rule featureCount_version:
    output: DATA+"/versions/subread.txt"
    conda: "envs/subread.yaml"
    shell:
        """
        featureCounts -v 2>&1 | tee {output}
        """    

rule make_strand_forward_txt:
    input:
        gtf=DATA+"/strand/targets.gtf",
        bam0=expand(DATA+"/bowtie2/{name}.bam",
                    name = SAMPLES_NAMES)[0]
    output: DATA+"/strand/forward.txt"
    conda: "envs/subread.yaml"
    shell:
        """
        featureCounts \
        -O -p --countReadPairs -B -P -C \
        -t gene -g gene_id -f \
        -s 1 \
        -a {input.gtf} \
        -o {output} \
	{input.bam0}
        """
        
rule make_strand_reverse_txt:
    input:
        gtf=DATA+"/strand/targets.gtf",
        bam0=expand(DATA+"/bowtie2/{name}.bam",
                    name = SAMPLES_NAMES)[0]
    output: DATA+"/strand/reverse.txt"
    conda: "envs/subread.yaml"
    shell:
        """
        featureCounts \
        -O -p --countReadPairs -B -P -C \
        -t gene -g gene_id -f \
        -s 2 \
        -a {input.gtf} \
        -o {output} \
	{input.bam0}
        """

if "orientation" in config:
    rule make_strand_sh:
        input:
        output: DATA+"/strand/results.sh"
        shell: "echo ORIENTATION="+config['orientation']+" >> {output}"
else:
    rule make_strand_sh:
        input:
            f=DATA+"/strand/forward.txt",
            r=DATA+"/strand/reverse.txt"
        output: DATA+"/strand/results.sh"
        shell:
            """
            {PIPELINE}/scripts/rnaseq-strand-analysis.pl \
            -o {output} \
            {input.f} \
            {input.r}
            """
# ------------------------------------------------------------------------
# Make Sinister and Naive profiles
# ------------------------------------------------------------------------

rule make_profiles:
    input: PROFILE_RESULTS

rule make_profile:
    input:
        bam=DATA+"/bowtie2/{name}.bam",
        orientation=DATA+"/strand/results.sh"
    output: DATA+"/profiles/.done.{name}"
    params:
        name="{name}",
    shell:
        """
        . {input.orientation}

        case X"$ORIENTATION"X in
            XforwardX) opt_r= ;;
            XreverseX) opt_r=-r ;;
            X*X) echo 1>&2 cannot happen ; exit 1
        esac

	samtools view -h {input.bam} \
	| {PIPELINE}/scripts/sam2profiles -2 $opt_r \
				 -e -s -n \
				 -d $(dirname {output}) \
				 -t {params.name}
        touch {output}
        """

# ------------------------------------------------------------------------
# Make SOUS profiles
# ------------------------------------------------------------------------

rule make_sous:
    input: DATA+"/inputs/genome.fna"
    output: DATA+"/sous/.done"
    params:
        outdir=DATA+"/sous"
    shell:
        """
        sous -d {params.outdir} -t sous -u uniqueuess {input}
        touch {output}
        """

# ------------------------------------------------------------------------
# Make count files
# ------------------------------------------------------------------------

rule make_counts:
    input: expand(DATA+"/counts/counts_{name}.txt", \
                  name = SAMPLES_NAMES)

rule make_counts_annotations_gtf:
    input: DATA+"/strand/targets.gtf",
    output: DATA+"/counts/annotations.gtf",
    shell: "cp -a {input} {output}"
    
rule run_featureCounts:
    input:
        gtf=DATA+"/counts/annotations.gtf",
        bam=DATA+"/bowtie2/{name}.bam",
        orientation=DATA+"/strand/results.sh"
    output: DATA+"/counts/counts_{name}.txt"
    conda: "envs/subread.yaml"
    shell:
        """
        . {input.orientation}

        case X"$ORIENTATION"X in
            XforwardX) s=1 ;;
            XreverseX) s=2 ;;
            X*X) echo 1>&2 cannot happen ; exit 1
        esac

        featureCounts \
        -O -p --countReadPairs -B -P -C -s $s \
        -t gene -g gene_id -f \
        -s 1 \
        -a {input.gtf} \
        -o {output} \
	{input.bam}
        """
        
# ------------------------------------------------------------------------
# Run DESeq2
# ------------------------------------------------------------------------

rule run_deseq2:
    input: DESEQ2_RESULTS

rule deseq2_version:
    output: DATA+"/versions/deseq2.txt"
    conda: "envs/deseq2.yaml"
    shell:
        """
cat <<EOF | Rscript - 2>&1 | tee {output}
suppressPackageStartupMessages(library("DESeq2"))
sessionInfo()
EOF
        """    

rule make_deseq2_regions_gtf:
    input: DATA+"/counts/annotations.gtf"
    output: DATA+"/deseq2/temp/regions.gtf"
    shell: "cp -a {input} {output}"


rule make_counts_txt:
    input:
        gtf=DATA+"/deseq2/temp/regions.gtf",
        counts=expand(DATA+"/counts/counts_{name}.txt", name=config['samples'].keys())
    output: DATA+"/deseq2/temp/counts_{tag}.txt"
    params:
        opt_t="{tag}",
        args = list(map(lambda name: str(name)+":"+DATA+"/counts/counts_"+str(name)+".txt",
                     config['samples'].keys()))
    shell:
        """
        {PIPELINE}/scripts/make-counts-table-from-featurecounts \
           {input.gtf} {params.args} > {output}
        """

rule make_deseq2_params_R:
    input: DATA+"/deseq2/temp/counts_{tag}.txt"
    output:
        params_r=DATA+"/deseq2/temp/params_{tag}.R",
        data_txt=DATA+"/deseq2/temp/data_{tag}.txt"
    params:
        opt_f=lambda wildcards: "-f "+str(config['exps'][wildcards.tag]["fdr"]) if "fdr" in config['exps'][wildcards.tag] else "",
        opt_p=lambda wildcards: "-p "+str(config['exps'][wildcards.tag]["pvalue"]) if "pvalue" in config['exps'][wildcards.tag] else "",
        opt_t=lambda wildcards: wildcards.tag,
        control_args=lambda wildcards: list(map((lambda name:config['exps'][wildcards.tag]['control_name']+":"+name),config['exps'][wildcards.tag]['control_samples'])),
        treatment_args=lambda wildcards: list(map((lambda name:config['exps'][wildcards.tag]['treatment_name']+":"+name),config['exps'][wildcards.tag]['treatment_samples'])),
    shell:
        """
        {PIPELINE}/scripts/prep-deseq2 -x -s {PIPELINE}/scripts \
	       -F parametric \
	       {params.opt_f} {params.opt_p} \
	       -d $(dirname {output.params_r}) \
	       -t _{params.opt_t} \
	       -c {input} \
	       {params.control_args} {params.treatment_args}
        """


rule run_deseq2_once:
    input:
        params_r=DATA+"/deseq2/temp/params_{tag}.R",
        data_txt=DATA+"/deseq2/temp/data_{tag}.txt"
    output:
        output_txt=DATA+"/deseq2/temp/output_{tag}.txt",
        xoutput_txt=DATA+"/deseq2/temp/output-extended_{tag}.txt",
    conda: "envs/deseq2.yaml"
    shell:
        """
        Rscript {PIPELINE}/scripts/run-deseq2 {input.params_r}
        """

rule make_results_txt:
    input:
        txt=DATA+"/deseq2/temp/output-extended_{tag}.txt",
        gtf=DATA+"/deseq2/temp/regions.gtf",        
    output: DATA+"/deseq2/results_{tag}.txt"
    params:
        aliases_txt = config['aliases_txt'] if 'aliases_txt' in config else []
    shell:
        """
        cat {input.txt} \
	| {PIPELINE}/scripts/deseq-output2 -t \
		     {input.gtf} \
		     {params.aliases_txt} \
		     > {output}
        """

rule make_results_gff:
    input:
        txt=DATA+"/deseq2/temp/output_{tag}.txt",
        gtf=DATA+"/deseq2/temp/regions.gtf",        
    output:
        DATA+"/deseq2/results_{tag}.gff"
    params:
        opt_c=lambda wildcards: "-c "+str(config['exps'][wildcards.tag]["foldchange"]) if "foldchange" in config['exps'][wildcards.tag] else "",
        aliases_txt = config['aliases_txt'] if 'aliases_txt' in config else []
    shell:
        """
        cat {input.txt} \
        | {PIPELINE}/scripts/deseq-output2 -g {params.opt_c} \
        	     {input.gtf} \
		     {params.aliases_txt} \
		     > {output}
        """


rule make_changed_gff:
    input: DATA+"/deseq2/results_{tag}.gff"
    output: DATA+"/deseq2/changed_{tag}.gff"
    shell:
        """
        set +e # if egrep matches nothing
        cat {input} | egrep '; colour [23];' > {output}
        set -e
        """


# ------------------------------------------------------------------------
# Get package versions
# ------------------------------------------------------------------------

rule make_versions:
    input: VERSIONS_RESULTS

# ------------------------------------------------------------------------
# Check Git status
# ------------------------------------------------------------------------

rule run_git:
    output: GIT_RESULTS
    shell:
        """
	(
	    cd {PIPELINE}
	    echo
	    ( set -x ; git status )
	    echo
	    ( set -x ; git log -n1 )
	) 2>&1 | tee {output}
        """ 
