from os.path import join


def list_samples():
    for sample in config["samples"].keys():
        yield sample


def get_from_sample(sample, what):
    out = config["samples"][sample][what]
    return out


# Directories & config
IONICE = config["ionice"]
RESULTS = config["results"]


# Wildcards
samples = list_samples()


rule all:
    input:
        expand(join(RESULTS, "output.model.txt"), sample=samples)


rule find_excluded_regions:
    """
    Regions to exclude from shuffling
    """
    input:
        peaks = lambda wildcards: get_from_sample(wildcards.sample, "peaks"),
    output:
        join(RESULTS, "regions_to_exclude.bed"),
    params:
        bl = config["bl_regions"],
    shell:
        """
        {IONICE} zless {input} {params} | cut -f 1-3 | sortBed | mergeBed \
            > {output}
        """

rule get_negative_bed:
    """
    List negative regions for GKM-SVM
    """
    input:
        peaks = lambda wildcards: get_from_sample(wildcards.sample, "peaks"),
        excl = rules.find_excluded_regions.output,
    output:
        join(RESULTS, "negative_regions.bed"),
    params:
        g = config["genome_size"],
        opts = "-i -noOverlapping -f 0 -seed 8789"
    shell:
        """
        {IONICE} bedtools shuffle {params.opts} -i {input.peaks} \
            -g {params.g} -excl {input.excl} > {output}
        """

rule get_fastas:
    """
    Make positive and negative files for GKM-SVM
    """
    input:
        peaks = lambda wildcards: get_from_sample(wildcards.sample, "peaks"),
        neg = rules.get_negative_bed.output,
    output:
        peaks = join(RESULTS, "peaks.fa"),
        neg = join(RESULTS, "negative_regions.fa"),
    params:
        fasta = config["genome"],
    shell:
        """
        {IONICE} bedtools getfasta -fi {params.fasta} -bed {input.peaks} \
            > {output.peaks}
        {IONICE} bedtools getfasta -fi {params.fasta} -bed {input.neg} \
            > {output.neg}
        """

rule train_gkm_svm:
    input:
        pos = rules.get_fastas.output.peaks,
        neg = rules.get_fastas.output.neg,
    output:
        join(RESULTS, "output.model.txt"),
    params:
        bin = join(config["gmkbin"], "gkmtrain"),
        handle = "output",
        opts = "-m {}".format(config["memory"])
    threads:
        config["cores"]
    shell:
        """
        /usr/bin/time -v {params.bin} -T {threads} {params.opts} \
            {input.pos} {input.neg} {params.handle}
        """