#!/usr/bin/env python3

import os
import pathlib2
import tempfile
from Bio import SeqIO
import multiprocessing

#############
# FUNCTIONS #
#############


def assembly_catalog_resolver(wildcards):
    if wildcards.name in assembly_catalog:
        return({'fasta': assembly_catalog[wildcards.name]})
    elif wildcards.name in polished_assemblies:
        return({'fasta': polished_assemblies[wildcards.name]})
    elif wildcards.name in merged_assemblies:
        return({'fasta': merged_assemblies[wildcards.name]})
    elif wildcards.name in final_assemblies:
        return({'fasta': final_assemblies[wildcards.name]})
    else:
        raise ValueError('missing {} in catalog'.format(wildcards.name))


def filter_fasta_by_length(input_fasta,
                           output_fasta,
                           length=10000):
    '''
    Read input_fasta file and write contigs longer than length to output_fasta
    '''
    SeqIO.write(sequences=(rec for rec in SeqIO.parse(input_fasta, 'fasta')
                           if len(rec) >= length),
                handle=output_fasta,
                format='fasta')


def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


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

###########
# GLOBALS #
###########

r1_raw = ['data/illumina/pe100/ASW_1.fastq.gz',
          'data/illumina/pe150/ASW_1.fastq.gz']
r2_raw = ['data/illumina/pe100/ASW_2.fastq.gz',
          'data/illumina/pe150/ASW_2.fastq.gz']
ont_raw = 'data/guppy_237.fastq.gz'
ont_tmp = 'output/000_tmp/guppy_237.fq'
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'
meraculous_config_file = 'src/meraculous_config.txt'
meraculous_threads = 72
dpon_ref = 'data/GCF_000355655.1_DendPond_male_1.0_genomic.fna'

# containers
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
flye_container = 'shub://TomHarrop/singularity-containers:flye_2.4'
canu_container = 'shub://TomHarrop/singularity-containers:canu_1.8'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
racon_container = 'shub://TomHarrop/singularity-containers:racon_1.3.2'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
quast_container = 'shub://TomHarrop/singularity-containers:quast_5.0.2'


# assembly catalog
assembly_catalog = {
    # 'flye_denovo': 'output/030_flye/de_novo/scaffolds.fasta',
    'meraculous': ('output/020_meraculous/k71_diplo2/'
                   'meraculous_final_results/final.scaffolds.fa'),
    'flye_denovo_full': 'output/030_flye/denovo_full/scaffolds.fasta',
    # 'flye_full_meraculous': ('output/040_merged_assemblies/'
    #                          'flye_full_meraculous/scaffolds.fasta'),
    'canu': 'output/035_canu/canu.contigs.fasta'
}

polished_assemblies = {
    'flye_denovo_full_polished': ('output/045_long_read_polishing/'
                                  'flye_denovo_full/'
                                  'flye_denovo_full.racon.fasta'),
    'flye_denovo_full_polished2': ('output/045_short_read_polishing/'
                                   'flye_denovo_full/'
                                   'flye_denovo_full.racon.fasta'),
    'meraculous_polished': ('output/045_long_read_polishing/meraculous/'
                            'meraculous.racon.fasta'),
    # 'meraculous_polished2': ('output/045_short_read_polishing/meraculous/'
    #                          'meraculous.racon.fasta'),
    'canu_polished': 'output/045_long_read_polishing/canu/canu.racon.fasta',
    'canu_polished2': 'output/045_short_read_polishing/canu/canu.racon.fasta'}

merged_assemblies = {
    # 'canu_flye':
    #     'output/040_merged_assemblies/canu_flye/scaffolds.fasta',
    # 'canu_flye_polished':
    #     'output/045_long_read_polishing/canu_flye/canu_flye.racon.fasta',
    # 'canu_flye_polished2':
    #     'output/045_short_read_polishing/canu_flye/canu_flye.racon.fasta'
}

final_assemblies = {
    'flye_denovo_full_final': ('output/047_final_polish/flye_denovo_full/'
                               'flye_denovo_full.racon.fasta'),
    'canu_final': ('output/047_final_polish/canu/'
                   'canu.racon.fasta')}


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
        expand('output/050_busco/run_{name}/full_table_{name}.tsv',
               name=list(assembly_catalog.keys()) +
               list(polished_assemblies.keys()) +
               list(final_assemblies.keys()))
        # 'output/055_quast/report.txt',
        # expand('output/050_busco/run_{name}/full_table_{name}.tsv',
        #        name=list(final_assemblies.keys()))



# general filtering rule
rule filter:
    input:
        fa = '{fa_name}.{fasta}'
    output:
        fa = '{fa_name}_filtered.{fasta}'
    params:
        length = 10000
    run:
        filter_fasta_by_length(input.fa, output.fa, params.length)


rule quast:
    input:
        assemblies = list(assembly_catalog.values()) +
        list(polished_assemblies.values()) +
        list(merged_assemblies.values()),
        ref = dpon_ref,
        pe = 'output/000_tmp/pe_reads.fq',
        ont = ont_tmp
    output:
        'output/055_quast/report.txt'
    log:
        'output/logs/055_quast/quast.log'
    params:
        labels = ','.join(list(assembly_catalog.keys()) +
                          list(polished_assemblies.keys()) +
                          list(merged_assemblies.keys())),
        outdir = 'output/055_quast'
    threads:
        meraculous_threads
    singularity:
        quast_container
    shell:
        'quast.py '
        '-r {input.ref} '
        '--threads {threads} '
        '--labels {params.labels} '
        '--eukaryote '
        '--large '
        '--k-mer-stats '
        '--circos '
        '--glimmer '
        '--rna-finding '
        '--conserved-genes-finding '    # BUSCO
        '--fragmented '
        '--pe12 {input.pe} '
        '--nanopore {input.ont} '
        '--output-dir {params.outdir} '
        '{input.assemblies} '
        '&> {log}'

rule busco:
    input:
        unpack(assembly_catalog_resolver),
        lineage = 'data/busco/endopterygota_odb9'
    output:
        'output/050_busco/run_{name}/full_table_{name}.tsv'
    log:
        resolve_path('output/logs/050_busco/busco_{name}.log')
    params:
        wd = 'output/050_busco',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage),
        tmpdir = tempfile.mkdtemp()
    threads:
        multiprocessing.cpu_count()
    priority:
        1
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path {params.tmpdir} '
        '--in {params.fasta} '
        '--out {wildcards.name} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species tribolium2012 '
        '--mode genome '
        '&> {log}'


# final polish
rule final_polish:
    input:
        # unpack(assembly_catalog_resolver),
        fasta = 'output/045_short_read_polishing/{name}/{name}.racon.fasta',
        aln = 'output/047_final_polish/{name}/aln.sam',
        fq = 'output/000_tmp/pe_reads.fq'
    output:
        'output/047_final_polish/{name}/{name}.racon.fasta'
    log:
        'output/logs/047_final_polish/{name}_racon.log'
    threads:
        multiprocessing.cpu_count()
    priority:
        0
    singularity:
        racon_container
    shell:
        'racon '
        '-t {threads} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'


rule map_final_polish:
    input:
        fasta = 'output/045_short_read_polishing/{name}/{name}.racon.fasta',
        fq = 'output/000_tmp/pe_reads.fq'
    output:
        'output/047_final_polish/{name}/aln.sam'
    params:
        prefix = 'output/047_final_polish/{name}/index'
    log:
        'output/logs/047_final_polish/{name}_bwa-mem.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log} '
        '; '
        'bwa mem '
        '-t {threads} '
        '-p '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2>> {log}'



# 04 wacky genome combinations + polishing
rule canu_flye:
    input:
        subassemblies = [
            'output/045_short_read_polishing/canu/canu.racon_filtered.fasta',
            ('output/045_short_read_polishing/'
             'flye_denovo_full/'
             'flye_denovo_full.racon.fasta')]
    output:
        'output/040_merged_assemblies/canu_flye/scaffolds.fasta'
    params:
        outdir = 'output/040_merged_assemblies/canu_flye',
        size = '1200m'
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/040_merged_assemblies/canu_flye.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--subassemblies {input.subassemblies} '
        '--iterations 0 '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


rule polish_short_reads:
    input:
        # unpack(assembly_catalog_resolver),
        fasta = 'output/045_long_read_polishing/{name}/{name}.racon.fasta',
        aln = 'output/045_short_read_polishing/{name}/aln.sam',
        fq = 'output/000_tmp/pe_reads.fq'
    output:
        'output/045_short_read_polishing/{name}/{name}.racon.fasta'
    log:
        'output/logs/045_short_read_polishing/{name}_racon.log'
    threads:
        multiprocessing.cpu_count()
    priority:
        0
    singularity:
        racon_container
    shell:
        'racon '
        '-t {threads} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'


rule polish_long_reads:
    input:
        unpack(assembly_catalog_resolver),
        aln = 'output/045_long_read_polishing/{name}/aln.sam',
        fq = ont_tmp,
    output:
        'output/045_long_read_polishing/{name}/{name}.racon.fasta'
    log:
        'output/logs/045_long_read_polishing/{name}_racon.log'
    threads:
        multiprocessing.cpu_count()
    priority:
        0
    singularity:
        racon_container
    shell:
        'racon '
        '-t {threads} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'


rule map_short_reads:
    input:
        fasta = 'output/045_long_read_polishing/{name}/{name}.racon.fasta',
        fq = 'output/000_tmp/pe_reads.fq'
    output:
        'output/045_short_read_polishing/{name}/aln.sam'
    params:
        prefix = 'output/045_short_read_polishing/{name}/index'
    log:
        'output/logs/045_short_read_polishing/{name}_bwa-mem.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log} '
        '; '
        'bwa mem '
        '-t {threads} '
        '-p '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2>> {log}'

rule map_long_reads:
    input:
        unpack(assembly_catalog_resolver),
        fq = ont_tmp
    output:
        'output/045_long_read_polishing/{name}/aln.sam'
    log:
        'output/logs/045_long_read_polishing/{name}_minimap.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax '
        'map-ont '
        '{input.fasta} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


rule flye_full_meraculous:
    input:
        subassemblies = [('output/020_meraculous/k71_diplo2/'
                          'meraculous_final_results/'
                          'final.scaffolds_filtered.fa'),
                         'output/030_flye/denovo_full/scaffolds_filtered.fasta']
    output:
        'output/040_merged_assemblies/flye_full_meraculous/scaffolds.fasta'
    params:
        outdir = 'output/040_merged_assemblies/flye_full_meraculous',
        size = '800m'
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/040_merged_assemblies/flye_full_meraculous.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--subassemblies {input.subassemblies} '
        '--iterations 0 '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


# 035 canu (since flye didn't work great)
rule canu:
    input:
        fq = ont_tmp
    output:
        'output/035_canu/canu.contigs.fasta'
    params:
        outdir = 'output/035_canu',
        size = '800m',
        prefix = 'canu'
    threads:
        multiprocessing.cpu_count()
    priority:
        5
    log:
        'output/logs/035_canu/canu.log'
    singularity:
        canu_container
    shell:
        'canu '
        '-d {params.outdir} '
        '-p {params.prefix} '
        'genomeSize={params.size} '
        'corMhapSensitivity=high '
        'corMinCoverage=0 '
        'corOutCoverage=100 '
        '-nanopore-raw {input.fq} '
        '&> {log}'


# 03 flye
rule flye_denovo_full:
    input:
        fq = ont_tmp
    output:
        'output/030_flye/denovo_full/scaffolds.fasta'
    params:
        outdir = 'output/030_flye/denovo_full',
        size = '800m'
    threads:
        multiprocessing.cpu_count()
    priority:
        10
    log:
        'output/logs/030_flye/denovo_full.log'
    singularity:
        flye_container
    shell:
        'flye '
        '--iterations 1 '
        # '--resume '
        '--nano-raw {input.fq} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'

rule flye:
    input:
        fq = ont_tmp
    output:
        'output/030_flye/de_novo/scaffolds.fasta'
    params:
        outdir = 'output/030_flye/de_novo',
        size = '800m'
    threads:
        multiprocessing.cpu_count()
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
        min(50, multiprocessing.cpu_count())
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
        threads = min(50, multiprocessing.cpu_count()),
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
        multiprocessing.cpu_count()
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

# speed up long read stuff with tmp unzip folder
rule tmp_unzip:
    input:
        ont_raw
    output:
        ont_tmp
    threads:
        3
    singularity:
        pigz_container
    shell:
        'pigz -dc '
        '{input} '
        '>{output}'

rule tmp_unzip_short:
    input:
        'output/010_trim-decon/pe_reads.fq.gz'
    output:
        'output/000_tmp/pe_reads.fq'
    threads:
        3
    singularity:
        pigz_container
    shell:
        'pigz -dc '
        '{input} '
        '>{output}'
