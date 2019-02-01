#!/usr/bin/env python3

import os
import pathlib2


#############
# FUNCTIONS #
#############


def busco_wildcard_resolver(wildcards):
    name_to_fasta = {
        'flye_denovo': 'output/030_flye/de_novo/scaffolds.fasta',
        'meraculous_filtered': '',
        'meraculous': ('output/020_meraculous/k71_diplo2/'
                       'meraculous_final_results/final.scaffolds.fa')
    }
    return({'fasta': name_to_fasta[wildcards.name]})


def write_config_file(fastq,
                      k,
                      diplo_mode,
                      dmin,
                      threads,
                      config_string,
                      config_file):
    '''
    Accept fastq file, threads config string and output location and write
    config
    '''
    print(fastq, threads, config_string, config_file)
    my_fastq = str(pathlib2.Path(fastq).resolve())
    my_conf = config_string.format(my_fastq, k, diplo_mode, dmin, threads)
    with open(config_file, 'wt') as f:
        f.write(my_conf)
    return True


def resolve_path(x):
    return str(pathlib2.Path(x).resolve())

###########
# GLOBALS #
###########

r1_raw = ['data/illumina/pe100/ASW_1.fastq.gz',
          'data/illumina/pe150/ASW_1.fastq.gz']
r2_raw = ['data/illumina/pe100/ASW_2.fastq.gz',
          'data/illumina/pe150/ASW_2.fastq.gz']
ont_raw = 'data/nanopore/merged_sorted.fq.gz'
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'
meraculous_config_file = 'src/meraculous_config.txt'
meraculous_threads = 64

# containers
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
flye_container = 'shub://TomHarrop/singularity-containers:flye_2.4'


########
# MAIN #
########

# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())


#########
# RULES #
#########

rule target:
    input:
        ('output/020_meraculous/k71_diplo2/'
         'meraculous_final_results/final.scaffolds.fa'),
        'output/030_flye/de_novo/scaffolds.fasta'

# 05 busco
rule busco_jobs:
    input:
        expand('output/050_busco/run_{name}/full_table_{name}.tsv',
               name=['flye_denovo',
                     'meraculous'])

rule busco:
    input:
        unpack(busco_wildcard_resolver),
        lineage = 'data/busco/endopterygota_odb9'
    output:
        'output/050_busco/run_{name}/full_table_{name}.tsv'
    log:
        resolve_path('output/logs/060_busco/busco_{name}.log')
    params:
        wd = 'output/050_busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage)
    threads:
        meraculous_threads
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path \"$(mktemp)\" '
        '--in {params.fasta} '
        '--out {wildcards.name} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species tribolium2012 '
        '--mode genome '
        '&> {log}'

# 04 wacky genome combinations + polishing


# 03 flye
rule flye:
    input:
        fq = ont_raw
    output:
        'output/030_flye/de_novo/scaffolds.fasta'
    params:
        outdir = 'output/030_flye/de_novo',
        size = '800m'
    threads:
        meraculous_threads
    log:
        'output/logs/030_flye/de_novo.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--asm-coverage 30 '
        '--nano-raw {input.fq} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'

# 02 meraculous
# rule filter_meraculous:


rule meraculous:
    input:
        fq = 'output/010_trim-decon/pe_reads.fq.gz',
        config = ('output/020_meraculous/'
                  'k{k}_diplo{diplo}/config.txt')
    output:
        ('output/020_meraculous/k{k}_diplo{diplo}/'
         'meraculous_final_results/final.scaffolds.fa')
    params:
        outdir = 'output/020_meraculous/k{k}_diplo{diplo}',
        dmin = '0'
    threads:
        min(50, meraculous_threads)
    log:
        'output/logs/020_meraculous/k{k}_diplo{diplo}.log'
    singularity:
        mer_container
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.config} '
        '-cleanup_level 2 '
        '&> {log}'

rule meraculous_config:
    input:
        fq = 'output/010_trim-decon/pe_reads.fq.gz'
    output:
        config = ('output/020_meraculous/'
                  'k{k}_diplo{diplo}/config.txt'),
    threads:
        1
    params:
        threads = min(50, meraculous_threads),
        dmin = '0'
    run:
        write_config_file(
            input.fq,
            wildcards.k,
            wildcards.diplo,
            params.dmin,
            params.threads,
            meraculous_config_string,
            output.config)

# 01 trim and decontaminate reads
rule trim_decon:
    input:
        r1 = 'output/010_trim-decon/r1.fq.gz',
        r2 = 'output/010_trim-decon/r2.fq.gz'
    output:
        fq = 'output/010_trim-decon/pe_reads.fq.gz',
        f_stats = 'output/010_trim-decon/filter-stats.txt',
        t_stats = 'output/010_trim-decon/trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/filter.log',
        trim = 'output/logs/010_trim-decon/trim.log',
        repair1 = 'output/logs/010_trim-decon/repair1.log',
        repair2 = 'output/logs/010_trim-decon/repair2.log'
    params:
        filter = bbduk_ref,
        trim = bbduk_adaptors
    threads:
        meraculous_threads // 2
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '2> {log.repair1} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'repair.sh '
        'in=stdin.fastq '
        'out={output.fq} '
        '2> {log.repair2} '

rule join_reads:
    input:
        r1 = r1_raw,
        r2 = r2_raw
    output:
        r1 = temp('output/010_trim-decon/r1.fq.gz'),
        r2 = temp('output/010_trim-decon/r2.fq.gz')
    singularity:
        bbduk_container
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'
