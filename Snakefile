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
flye = 'shub://TomHarrop/assemblers:flye_2.6'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'

#########
# RULES #
#########

rule target:
    input:
        'output/020_flye-assemble/assembly.fasta'

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
rule flye_assemble:
    input:
        pool = 'output/010_raw/pool.fq',
        asw47 = 'output/010_raw/asw47.fq'
    output:
        'output/020_flye-assemble/assembly.fasta'
    params:
        outdir = 'output/020_flye_assemble',
        size = '1.2g'
    threads:
        max(128, multiprocessing.cpu_count())
    log:
        'output/logs/020_flye-assemble.log'
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

