from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

snv_indel_vcf = sys.argv[1]
sv_vcf = sys.argv[2]
ped_uri = sys.argv[3]
prefix = sys.argv[4]
omim_uri = sys.argv[5]
sv_gene_fields = sys.argv[6].split(',')  # ['PREDICTED_LOF', 'PREDICTED_INTRAGENIC_EXON_DUP']
build = sys.argv[7]
cores = sys.argv[8]  # string
mem = int(np.floor(float(sys.argv[9])))

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": "2",
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
        #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
        #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def filter_mt(mt):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                    'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
    mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
        hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter by Impact and splice/noncoding consequence
    splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
    keep_vars = ['non_coding_transcript_exon_variant']
    mt = mt.filter_rows(
        (hl.set(splice_vars + keep_vars).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
        (hl.array(['HIGH', 'MODERATE']).contains(
        mt.vep.transcript_consequences.IMPACT))
        )
    return mt 

## STEP 1: Merge SNV/Indel VCF with SV VCF (or just one of them)
# Load SNV/Indel VCF
if snv_indel_vcf!='NA':
    snv_mt = hl.import_vcf(snv_indel_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)
    csq_columns = hl.get_vcf_metadata(snv_indel_vcf)['info']['CSQ']['Description'].split('Format: ')[1].split('|')

    snv_mt = snv_mt.annotate_rows(vep=snv_mt.info)
    transcript_consequences = snv_mt.vep.CSQ.map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                        {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                            for i, col in enumerate(csq_columns)})))

    snv_mt = snv_mt.annotate_rows(vep=snv_mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
    snv_mt = snv_mt.annotate_rows(vep=snv_mt.vep.select('transcript_consequences'))

    # filter SNV/Indel MT
    snv_mt = snv_mt.explode_rows(snv_mt.vep.transcript_consequences)
    snv_mt = filter_mt(snv_mt)

    snv_gene_field_str = 'vep.transcript_consequences.SYMBOL'
    snv_gene_field = None

    for i, field in enumerate(snv_gene_field_str.split('.')):
        if i==0:
            snv_gene_field = snv_mt[field]
        else:
            snv_gene_field = snv_gene_field[field]
    snv_mt = snv_mt.annotate_rows(gene=snv_gene_field)
    snv_mt = snv_mt.annotate_rows(variant_type='SNV/Indel')

# Load SV VCF
if sv_vcf!='NA':
    sv_mt = hl.import_vcf(sv_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

    sv_mt = sv_mt.annotate_rows(gene=hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in sv_gene_fields]))))
    sv_mt = sv_mt.annotate_rows(variant_type='SV')

    sv_mt = sv_mt.explode_rows(sv_mt.gene)

if (snv_indel_vcf!='NA') and (sv_vcf!='NA'):
    sv_info_fields, sv_entry_fields = list(sv_mt.row.info), list(sv_mt.entry)
    snv_info_fields, snv_entry_fields = list(snv_mt.row.info), list(snv_mt.entry)

    sv_missing_entry_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_entry_fields, sv_entry_fields)}
    snv_missing_entry_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_entry_fields, snv_entry_fields)}

    sv_missing_info_fields = {field: str(snv_mt.info[field].dtype) for field in np.setdiff1d(snv_info_fields, sv_info_fields)}
    snv_missing_info_fields = {field: str(sv_mt.info[field].dtype) for field in np.setdiff1d(sv_info_fields, snv_info_fields)}

    sv_mt = sv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in sv_missing_entry_fields.items()})
    snv_mt = snv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in snv_missing_entry_fields.items()})

    sv_mt = sv_mt.select_entries(*sorted(list(sv_mt.entry)))
    snv_mt = snv_mt.select_entries(*sorted(list(snv_mt.entry)))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in sv_missing_info_fields.items()}))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in snv_missing_info_fields.items()}))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.select(*sorted(list(sv_mt.info))))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.select(*sorted(list(snv_mt.info))))

    # VEP
    snv_vep_fields = {field: str(snv_mt.vep.transcript_consequences[field].dtype) for field in list(snv_mt.row.vep.transcript_consequences)}
    sv_mt = sv_mt.annotate_rows(vep=hl.struct(transcript_consequences=
            {field: hl.missing(dtype) for field, dtype in snv_vep_fields.items()}))

    # Annotate OMIM in SVs
    omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
    sv_mt = sv_mt.key_rows_by('gene')
    sv_mt = sv_mt.annotate_rows(vep=sv_mt.vep.annotate(
        transcript_consequences=sv_mt.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].inheritance_code, ''))))

    sv_mt = sv_mt.key_rows_by().select_rows(*sorted(list(sv_mt.row))).key_rows_by('locus','alleles')
    snv_mt = snv_mt.key_rows_by().select_rows(*sorted(list(snv_mt.row))).key_rows_by('locus','alleles')

    # Subset shared samples 
    sv_samps = sv_mt.s.collect()
    snv_samps = snv_mt.s.collect()
    shared_samps = list(np.intersect1d(sv_samps, snv_samps))

    sv_mt = sv_mt.filter_cols(hl.array(shared_samps).contains(sv_mt.s))
    snv_mt = snv_mt.filter_cols(hl.array(shared_samps).contains(snv_mt.s))

    variant_types = 'SV_SNV_Indel'
    merged_mt = sv_mt.key_cols_by().union_rows(snv_mt.key_cols_by()).key_cols_by('s')

elif snv_indel_vcf!='NA':
    variant_types = 'SNV_Indel'
    merged_mt = snv_mt

elif sv_vcf!='NA':
    variant_types = 'SV'
    merged_mt = sv_mt


## STEP 2: Get CompHets
pedigree = hl.Pedigree.read(ped_uri, delimiter='\t')
trio_samples = list(np.array([[trio.s, trio.pat_id, trio.mat_id] for trio in pedigree.complete_trios()]).flatten())

# Mendel errors
def get_mendel_errors(mt, phased_tm):
    all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
    phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)
    return phased_tm

def phase_by_transmission_aggregate_by_gene(tm, mt):
    # filter out calls that are hom ref in proband
    tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())

    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')
    phased_tm = get_mendel_errors(mt, phased_tm)
    gene_phased_tm = phased_tm.key_rows_by('locus','alleles','gene')

    gene_agg_phased_tm = (gene_phased_tm.group_rows_by(gene_phased_tm.gene)
        .aggregate_rows(locus_alleles = hl.agg.collect(gene_phased_tm.row_key),
                       variant_type = hl.agg.collect(gene_phased_tm.variant_type))
        .aggregate_entries(proband_PBT_GT = hl.agg.collect(gene_phased_tm.proband_entry.PBT_GT).filter(hl.is_defined),
                          proband_GT = hl.agg.collect(gene_phased_tm.proband_entry.GT).filter(hl.is_defined))).result()
    return gene_phased_tm, gene_agg_phased_tm

def get_subset_tm(mt, samples, keep=True, complete_trios=False):
    subset_mt = mt.filter_cols(hl.array(samples).contains(mt.s), keep=keep)

    # remove variants missing in subset samples
    subset_mt = hl.variant_qc(subset_mt)
    subset_mt = subset_mt.filter_rows(subset_mt.variant_qc.AC[1]>0)
    subset_mt = subset_mt.drop('variant_qc')

    subset_tm = hl.trio_matrix(subset_mt, pedigree, complete_trios=complete_trios)
    return subset_tm

def get_non_trio_comphets(mt):
    non_trio_tm = get_subset_tm(mt, trio_samples, keep=False)
    non_trio_gene_phased_tm, non_trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(non_trio_tm, mt)

    # different criteria for non-trios
    potential_comp_hets_non_trios = non_trio_gene_agg_phased_tm.filter_rows(
            hl.agg.count_where(non_trio_gene_agg_phased_tm.proband_GT.size()>1)>0
    )
          
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.explode_rows(potential_comp_hets_non_trios.locus_alleles)
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.key_rows_by(potential_comp_hets_non_trios.locus_alleles.locus, potential_comp_hets_non_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_non_trios = potential_comp_hets_non_trios.filter_entries(potential_comp_hets_non_trios.proband_GT.size()>1)
    
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.key_rows_by('locus', 'alleles', 'gene')
    non_trio_gene_phased_tm = non_trio_gene_phased_tm.annotate_entries(proband_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_GT),
                                                                      proband_PBT_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_gene_phased_tm.row_key, non_trio_gene_phased_tm.col_key].proband_PBT_GT))

    gene_phased_tm_comp_hets_non_trios = non_trio_gene_phased_tm.semi_join_rows(potential_comp_hets_non_trios.rows()).key_rows_by('locus', 'alleles')
    return gene_phased_tm_comp_hets_non_trios

def get_trio_comphets(mt):
    trio_tm = get_subset_tm(mt, trio_samples, keep=True, complete_trios=True)

    trio_gene_phased_tm, trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(trio_tm, mt)

    # different criteria for trios (requires phasing)
    potential_comp_hets_trios = trio_gene_agg_phased_tm.filter_rows(
        hl.agg.count_where(hl.set(trio_gene_agg_phased_tm.proband_PBT_GT).size()>1)>0
    )
    potential_comp_hets_trios = potential_comp_hets_trios.explode_rows(potential_comp_hets_trios.locus_alleles)
    potential_comp_hets_trios = potential_comp_hets_trios.key_rows_by(potential_comp_hets_trios.locus_alleles.locus, potential_comp_hets_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_trios = potential_comp_hets_trios.filter_entries(hl.set(potential_comp_hets_trios.proband_PBT_GT).size()>1)

    trio_gene_phased_tm = trio_gene_phased_tm.key_rows_by('locus', 'alleles', 'gene')
    trio_gene_phased_tm = trio_gene_phased_tm.annotate_entries(proband_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_GT),
                                                               proband_PBT_GT_set=hl.set(
        potential_comp_hets_trios[trio_gene_phased_tm.row_key, trio_gene_phased_tm.col_key].proband_PBT_GT))

    trio_gene_phased_tm = trio_gene_phased_tm.filter_entries(trio_gene_phased_tm.proband_PBT_GT_set.size()>1)  # this actually seems necessary but idk why tbh
    gene_phased_tm_comp_hets_trios = trio_gene_phased_tm.semi_join_rows(potential_comp_hets_trios.rows()).key_rows_by('locus', 'alleles')
    return gene_phased_tm_comp_hets_trios

def get_transmission(phased_tm_ht):
    phased_tm_ht = phased_tm_ht.annotate(transmission=hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm_ht.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm_ht

# Get CompHets
merged_trio_comphets = get_trio_comphets(merged_mt)
merged_non_trio_comphets = get_non_trio_comphets(merged_mt)
merged_comphets = merged_trio_comphets.entries().union(merged_non_trio_comphets.entries())


# XLR only
merged_tm = hl.trio_matrix(merged_mt, pedigree, complete_trios=False)
gene_phased_tm, gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(merged_tm, merged_mt)
xlr_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4'))
xlr_phased_tm = xlr_phased_tm.filter_entries((xlr_phased_tm.proband_entry.GT.is_non_ref()) &
                            (~xlr_phased_tm.is_female)).key_rows_by('locus', 'alleles')

# HomVar in proband only
phased_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
phased_hom_var = phased_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(phased_hom_var.proband_entry.GT))>0).key_rows_by('locus', 'alleles')

xlr_phased_tm = xlr_phased_tm.annotate_rows(variant_type='XLR')
phased_hom_var = phased_hom_var.annotate_rows(variant_type='hom_var')
merged_comphets = merged_comphets.annotate(variant_type='comphet')

merged_comphets_xlr_hom_var = merged_comphets.drop('proband_GT_set','proband_PBT_GT_set').union(xlr_phased_tm.entries()).union(phased_hom_var.entries())
merged_comphets_xlr_hom_var = get_transmission(merged_comphets_xlr_hom_var)

merged_comphets_xlr_hom_var.entries().flatten().export(f"{prefix}_{variant_types}_comp_hets_xlr_hom_var.tsv.gz")