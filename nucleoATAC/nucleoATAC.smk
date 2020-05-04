from os.path import join


def list_samples:
    for sample in config["sample"].keys():
        yield sample


def get_from_sample(sample, what):
    out = config["samples"][sample][what]
    return out


RESULTS = config["results"]


rule all:
    expand(join(RESULTS, "{sample}.nucpos.bed.gz"), sample=list_samples())


rule run:
    input:
        bam = lambda wildcards: get_from_sample(wildcards.sample, "bam"),
        peaks = lambda wildcards: get_from_sample(wildcards.sample, "peaks"),
    output:
        join(RESULTS, "{sample}.nucpos.bed.gz"),
    params:
        genome = config["genome"],
        handle = "{sample}",
    threads:
        config["cores"]
    shell:
        """
        nucleoatac run --bed {input.peaks} --bam {input.bam} \
            --fasta {params.genome} --out {params.handle} --cores {threads}
        """