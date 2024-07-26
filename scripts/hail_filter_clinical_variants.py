from pyspark.sql import SparkSession
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

vcf_file = sys.argv[1]
prefix = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))
ped_uri = sys.argv[5]
ac_threshold = int(sys.argv[6])
gnomad_af_threshold = float(sys.argv[7])

def get_transmission(df):
    '''
    df: trio matrix (tm) phased with PBT_GT converted to Pandas DataFrame
    returns Pandas Series
    ''' 
    return df['proband_entry.PBT_GT'].astype(str).map({
        '0|0': 'uninherited', 
        '0|1': 'inherited_from_mother', 
        '1|0': 'inherited_from_father', 
        '1|1': 'inherited_from_both', 
        'None': 'unknown'
    })

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

mt = hl.import_vcf(vcf_file, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)

header = hl.get_vcf_metadata(vcf_file)
csq_columns = header['info']['CSQ']['Description'].split('Format: ')[1].split('|')

# split VEP CSQ string
mt = mt.annotate_rows(vep=mt.info)
transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                       {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                        for i, col in enumerate(csq_columns)}), 
                                                        hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                        for i, col in enumerate(csq_columns)})))

mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))

mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)),  
                             gnomADg_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0]!='', 
                                                  hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADg_AF))[0])),
                             gnomADe_AF=hl.or_missing(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0]!='', 
                                                  hl.float(hl.array(hl.set(mt.vep.transcript_consequences.gnomADe_AF))[0])))
mt = mt.annotate_rows(gnomad_af=hl.max([mt.gnomADg_AF, mt.gnomADe_AF]))

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Output 1: grab ClinVar only
clinvar_tm = phased_tm.filter_rows((phased_tm.info.CLNSIG[0].matches('Pathogenic') | phased_tm.info.CLNSIG[0].matches('pathogenic')))
clinvar_tm = clinvar_tm.filter_entries((clinvar_tm.proband_entry.GT.is_non_ref()) | 
                                   (clinvar_tm.mother_entry.GT.is_non_ref()) |
                                   (clinvar_tm.father_entry.GT.is_non_ref()))
clinvar_tm = clinvar_tm.annotate_rows(variant_type='ClinVar_P/LP')
# clinvar_df = clinvar_tm.entries().to_pandas()
# clinvar_df['transmission'] = get_transmission(clinvar_df)

# filter out ClinVar benign
mt = mt.filter_rows((hl.is_missing(mt.info.CLNSIG)) |
    ~(mt.info.CLNSIG[0].matches('Benign') | mt.info.CLNSIG[0].matches('benign')))

# filter PASS
mt = mt.filter_rows(mt.filters.size()==0)

# filter out variants containing only these consequences
exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']

mt = mt.filter_rows(hl.set(exclude_csqs).intersection(mt.all_csqs).size()!=mt.all_csqs.size())

# filter by AC and gnomAD AF
mt = mt.filter_rows(mt.info.cohort_AC<=ac_threshold)
mt = mt.filter_rows((mt.gnomad_af<=gnomad_af_threshold) | (hl.is_missing(mt.gnomad_af)))

# export intermediate VCF
hl.export_vcf(mt, prefix+'_clinical.vcf.bgz', metadata=header)

# export ClinVar TSV
clinvar_tm.entries().export(prefix+'_clinvar_variants.tsv.gz', delimiter='\t')
# clinvar_df.to_csv(prefix+'_clinvar_variants.tsv.gz', sep='\t', index=False)
