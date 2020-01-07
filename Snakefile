#!/usr/bin/env python3

from pathlib import Path
import psutil
import tempfile

#############
# FUNCTIONS #
#############

def busco_target_resolver(wildcards):
    return {'fasta': busco_targets[wildcards.assembly]}


def resolve_path(path):
    return str(Path(path).resolve())


###########
# GLOBALS #
###########

pool_raw = 'data/ont-reads/pool.fq.gz'
asw47_raw = 'data/ont-reads/asw47.fq.gz'

# containers
busco = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
flye = 'shub://TomHarrop/assemblers:flye_2.6-g47548b8'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
purge_haplotigs = 'shub://TomHarrop/assembly-utils:purge_haplotigs_0b9afdf'

# resources
cpus = psutil.cpu_count()
mem = psutil.virtual_memory().total
fraction_to_use = 0.75
gb_mem_per_thread = round((mem * fraction_to_use) / cpus
                          / (1024 * 1024 * 1024))

# busco jobs
busco_targets = {
    'flye_polish': 'output/025_flye-polish/assembly.fasta',
    'flye_assemble': 'output/020_flye-assemble/assembly.fasta',
    'flye_asw47': 'output/027_flye-asw47/assembly.fasta'
}


#########
# RULES #
#########

rule target:
    input:
        expand('output/099_busco/run_{assembly}/full_table_{assembly}.tsv',
               assembly=list(busco_targets.keys())),
        'output/030_purge-haplotigs/histogram.png'

# busco
rule busco:
    input:
        unpack(busco_target_resolver),
        lineage = 'data/busco/endopterygota_odb9'
    output:
        'output/099_busco/run_{assembly}/full_table_{assembly}.tsv'
    log:
        resolve_path('output/logs/busco.{assembly}.log')
    params:
        wd = 'output/099_busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage),
        tmpdir = tempfile.mkdtemp()
    threads:
        multiprocessing.cpu_count() // 2
    singularity:
        busco
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path {params.tmpdir} '
        '--in {params.fasta} '
        '--out {wildcards.assembly} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species tribolium2012 '
        '--mode genome '
        '&> {log}'

# 03 purge haplotigs
rule haplotigs_hist:
    input:
        assembly = 'output/030_purge-haplotigs/ref.fasta',
        fai = 'output/030_purge-haplotigs/ref.fasta.fai',
        bam = 'output/030_purge-haplotigs/aligned.bam'
    output:
        'output/030_purge-haplotigs/histogram.png'
    params:
        wd = 'output/030_purge-haplotigs',
        assembly = lambda wildcards, input: resolve_path(input.assembly),
        bam = lambda wildcards, input: resolve_path(input.bam)
    log:
        resolve_path('output/logs/haplotigs_hist.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        purge_haplotigs
    shell:
        'cd {params.wd} || exit 1 ; '
        'purge_haplotigs '
        '-b {params.bam} '
        '-g {params.assembly} '
        '-t {threads} '
        '&> {log}'


rule haplotigs_sort:
    input:
        'output/030_purge-haplotigs/aligned.sam'
    output:
        'output/030_purge-haplotigs/aligned.bam'
    log:
        'output/logs/haplotigs_sort.log'
    singularity:
        purge_haplotigs
    shell:
        'samtools sort '
        '-m {gb_mem_per_thread}G '
        '-o {output} {input} '
        '-@ {cpus} ; '
        'samtools index {output}'

rule haplotigs_map:
    input:
        assembly = 'output/030_purge-haplotigs/ref.fasta',
        short_reads = 'data/short_reads.fastq'
    output:
        pipe = pipe('output/030_purge-haplotigs/aligned.sam')
    log:
        'output/logs/haplotigs_map.log'
    threads:
        multiprocessing.cpu_count() - 1
    singularity:
        purge_haplotigs
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax sr '
        '--secondary=no '
        '{input.assembly} '
        '{input.short_reads} '
        '>> {output.pipe} '
        '2> {log}'

rule haplotigs_ref:
    input:
        'output/025_flye-polish/assembly.fasta'
    output:
        fa = 'output/030_purge-haplotigs/ref.fasta',
        fai = 'output/030_purge-haplotigs/ref.fasta.fai'
    singularity:
        purge_haplotigs
    shell:
        'cp {input} {output.fa} '
        '; '
        'samtools faidx {output.fa}'

# 02 flye
rule flye_polish:
    input:
        assembly = 'output/020_flye-assemble/assembly.fasta',
        asw47 = 'output/010_raw/asw47.fq'
    output:
        'output/025_flye-polish/assembly.fasta'
    params:
        indir = 'output/020_flye-assemble',
        outdir = 'output/025_flye-polish',
        size = '1.2g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/flye_polish.log'
    singularity:
        flye
    shell:
        'cp -r {params.indir}/* {params.outdir} ; '
        'flye '
        '--nano-raw {input.asw47}  '
        '--resume-from polishing '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


rule flye_assemble:
    input:
        pool = 'output/010_raw/pool.fq',
        asw47 = 'output/010_raw/asw47.fq'
    output:
        'output/020_flye-assemble/assembly.fasta'
    params:
        outdir = 'output/020_flye-assemble',
        size = '1.2g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/flye_assemble.log'
    singularity:
        flye
    shell:
        'flye '
        '--nano-raw {input}  '
        '--iterations 0 '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'

rule flye_asw47:
    input:
        asw47 = 'output/010_raw/asw47.fq'
    output:
        'output/027_flye-asw47/assembly.fasta'
    params:
        outdir = 'output/027_flye-asw47',
        size = '1.2g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/flye_asw47.log'
    singularity:
        flye
    shell:
        'flye '
        '--nano-raw {input}  '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


# porechop
rule remove_ont_adaptors:
    input:
        'data/ont-reads/{group}.fq.gz'
    output:
        'output/010_raw/{group}.fq'
    log:
        'output/logs/remove_ont_adaptors.{group}.log'
    threads:
        multiprocessing.cpu_count() // 2
    singularity:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--check_reads 1000 '
        '--discard_middle '
        '&> {log}'

