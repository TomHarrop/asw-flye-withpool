#!/usr/bin/env python3

import multiprocessing

#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########

pool_raw = 'data/ont-reads/pool.fq.gz'
asw47_raw = 'data/ont-reads/asw47.fq.gz'

# containers
flye = 'shub://TomHarrop/assemblers:flye_2.6-g47548b8'
# flye = 'shub://TomHarrop/assemblers:flye_2.6'
# flye = 'shub://TomHarrop/singularity-containers:flye_2.5'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'

#########
# RULES #
#########

rule target:
    input:
        'output/025_flye-polish/assembly.fasta'

# rule busco:
#     input:
#         unpack(assembly_catalog_resolver),
#         lineage = 'data/busco/endopterygota_odb9'
#     output:
#         'output/050_busco/run_{name}/full_table_{name}.tsv'
#     log:
#         resolve_path('output/logs/050_busco/busco_{name}.log')
#     params:
#         wd = 'output/050_busco',
#         fasta = lambda wildcards, input: resolve_path(input.fasta),
#         lineage = lambda wildcards, input: resolve_path(input.lineage),
#         tmpdir = tempfile.mkdtemp()
#     threads:
#         multiprocessing.cpu_count()
#     priority:
#         1
#     singularity:
#         busco_container
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'run_BUSCO.py '
#         '--force '
#         '--tmp_path {params.tmpdir} '
#         '--in {params.fasta} '
#         '--out {wildcards.name} '
#         '--lineage {params.lineage} '
#         '--cpu {threads} '
#         '--species tribolium2012 '
#         '--mode genome '
#         '&> {log}'


# to polish this:
# https://github.com/fenderglass/Flye/issues/98


# 03 flye
rule flye_polish:
    input:
        assembly = 'output/025_flye-polish',
        asw47 = 'output/010_raw/asw47.fq'
    output:
        'output/025_flye-polish/assembly.fasta'
    params:
        outdir = 'output/025_flye-polish',
        size = '1.2g'
    threads:
        min(128, multiprocessing.cpu_count())
    log:
        'output/logs/flye_polish.log'
    singularity:
        flye
    shell:
        'flye '
        '--nano-raw {input.asw47}  '
        '--resume-from polishing '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


rule flye_duplicate_assembly:
    input:
        'output/020_flye-assemble/assembly.fasta'
    output:
        directory('output/025_flye-polish')
    params:
        wd = 'output/020_flye-assemble'
    singularity:
        flye
    shell:
        'cp -r {params.wd} {output}'


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

