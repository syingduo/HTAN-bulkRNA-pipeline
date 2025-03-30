import csv
from collections import namedtuple
from textwrap import dedent
from pathlib import Path
import pandas as pd
import logging


# The mapping of sample name to other information.
BAM_MAP_PTH = config['bam_map']
GENE_GTF_PTH = config['gdc_gtf']  # gencode.v22.annotation.gtf
GENE_INFO_PTH = config['gdc_gene_info']  # gencode.gene.info.v22.tsv
WORKFLOW_ROOT = config['workflow_root']  # Path to this repository
FC_STRAND = config['fc_strand']   # featureCounts strandness option
# Useful to change the BAM map in Docker
REPLACE_BAM_PTH = config.get('replace_bam_path', None)

_logger = logging.getLogger(__name__)

# Read all the cases to process
SAMPLES = set(pd.read_csv(BAM_MAP_PTH, sep='\t')['HTAN_Specimen_ID'])

# Select all the available samples of the selected cases.
bam_pth = {}
with open(BAM_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect="excel-tab")
    for row in reader:
        sample_id = row["HTAN_Specimen_ID"]
        bam_pth[sample_id] = Path(row["Path"])

def get_bam_inputs(wildcards):
    bam = bam_pth[wildcards.sample]
    return {
        "bam": bam,
        "bai": str(bam) + ".bai"
    }

rule featurecounts_readcount:
    """Readcount by featureCounts."""
    output: count_tsv=temp('featurecounts_readcount/{sample}.tsv')
    input: unpack(get_bam_inputs)
    log: 'logs/featurecounts/{sample}.log'
    params:
        gtf=GENE_GTF_PTH,
        strand=FC_STRAND
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 16000 + 16000 * (attempt - 1)
    threads: 16
    group: "featurecounts"
    shell:
        'featureCounts '
        '-g gene_id '  # feature id (-i in htseq)
        '-t exon '  # feature type (-t in htseq)
        '-T {threads} '
        '-Q 10 '  # htseq set this minimal mapping quality by default
        '-p '  # pair-end reads are considered one fragment; default HTSeq behavior
        '-B '  # both reads of a read pair need to be mapped
        '-a {params.gtf} '
        '-s {params.strand} '
        '-o {output.count_tsv} {input.bam} 2> {log}'


rule compress_featurecounts:
    """Shrink and compress featureCounts output."""
    output: 'featurecounts_readcount/{sample}.tsv.gz'
    input: rules.featurecounts_readcount.output.count_tsv
    threads: 2
    group: "featurecounts"
    shell: 'python {WORKFLOW_ROOT}/expression_quant/shrink_featurecounts.py {input} | gzip -9 -c > {output}'


rule generate_fpkm:
    """Generate FPKM, FPKM-UQ and TPM from the readcount."""
    output: fpkm='standardized_readcount/{sample}.tsv.gz'
    input: rc=rules.compress_featurecounts.output[0],
           gene_info=GENE_INFO_PTH
    shell: 'python {WORKFLOW_ROOT}/expression_quant/gen_fpkm.py {input.gene_info} {input.rc} {output.fpkm}'


rule all_featurecounts_readcount:
    input:
        counts=expand(rules.compress_featurecounts.output[0], sample=SAMPLES)


rule all_fpkms:
    input: fpkms=expand(rules.generate_fpkm.output.fpkm, sample=SAMPLES)
