version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterClinicalVariantsSV {
    input {
        File vcf_file
        File ped_uri
        File bed_files_with_header
        String cohort_prefix
        String genome_build='GRCh38'
        String hail_docker
        String variant_interpretation_docker

        Float bed_overlap_threshold=0.5
        Float gnomad_af_threshold=0.05
    }

    call vcfToBed {
        input:
        vcf_file=vcf_file,
        cohort_prefix=cohort_prefix,
        variant_interpretation_docker=variant_interpretation_docker
    }

    Array[Pair[String, File]] bed_files_map = read_map(bed_files_with_header)

    scatter (pair in bed_files_map) {
        call intersectBed {
            input:
            bed_file=vcfToBed.bed_output,
            ref_bed_with_header=pair.right,
            cohort_prefix=cohort_prefix,
            bed_overlap_threshold=bed_overlap_threshold,
            variant_interpretation_docker=variant_interpretation_docker
        }
    }

    call annotateVCFWithBeds {
        input:
        vcf_file=vcf_file,
        intersect_bed_files=intersectBed.intersect_bed,
        bed_files_with_header=bed_files_with_header,
        cohort_prefix=cohort_prefix,
        genome_build=genome_build,
        hail_docker=hail_docker
    }

    call filterVCF {
        input:
        vcf_file=annotateVCFWithBeds.annotated_vcf,
        ped_uri=ped_uri,
        genome_build=genome_build,
        hail_docker=hail_docker,
        gnomad_af_threshold=gnomad_af_threshold
    }

    output {
        File pathogenic_tsv = filterVCF.pathogenic_tsv
        File filtered_vcf = filterVCF.filtered_vcf
        File filtered_vcf_idx = filterVCF.filtered_vcf_idx
    }
}

task vcfToBed {
    input{
        File vcf_file
        String cohort_prefix
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        svtk vcf2bed --no-samples ~{vcf_file} ~{cohort_prefix}.bed.gz
    >>>

    output {
        File bed_output = "~{cohort_prefix}.bed.gz"
    }
}

task intersectBed {
    input {
        File bed_file
        File ref_bed_with_header
        String cohort_prefix
        String variant_interpretation_docker

        Float bed_overlap_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([bed_file, ref_bed_with_header], 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String ref_bed_with_header_str = basename(ref_bed_with_header, '.bed')

    command <<<
        set -eou pipefail
        tail -n +2 ~{ref_bed_with_header} > ref.bed 
        bedtools intersect -wao -f ~{bed_overlap_threshold} -r -a ~{bed_file} -b ref.bed | bgzip > ~{cohort_prefix}_~{ref_bed_with_header_str}.bed.gz
    >>>

    output {
        File intersect_bed = "~{cohort_prefix}_~{ref_bed_with_header_str}.bed.gz"
    }
}

task annotateVCFWithBeds {
    input {
        File vcf_file
        Array[File] intersect_bed_files
        File bed_files_with_header
        String cohort_prefix
        String genome_build
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, 'GB') + size(intersect_bed_files, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    set -eou pipefail
    cat <<EOF > annotate_vcf.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    intersect_bed_files = sys.argv[2].split(',')
    bed_files_with_header = sys.argv[3]
    cohort_prefix = sys.argv[4]
    genome_build = sys.argv[5]
    cores = sys.argv[6]
    mem = int(np.floor(float(sys.argv[7])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    header = hl.get_vcf_metadata(vcf_file)

    ref_bed_df = pd.read_csv(bed_files_with_header, sep='\t', header=None)
    annot_names = ref_bed_df[0].tolist()
    ref_bed_files = ref_bed_df[1].tolist()

    def annotate_with_bed(intersect_bed, annot_name, ref_bed_with_header_uri, mt, header):
        overlap_bed = hl.import_table(intersect_bed, force_bgz=True, no_header=True, types={f"f{i}": 'int' for i in [1,2,6,7]})
        overlap_bed = overlap_bed.annotate(f1=overlap_bed.f1 + 1,  # adjust for bed 0-based coordinates
                                        f6=overlap_bed.f6 + 1)

        fields = list(overlap_bed.row)
        overlap_field = fields[-1]

        overlap_bed = overlap_bed.annotate(sv_len=overlap_bed.f2-overlap_bed.f1, 
                            ref_len=overlap_bed.f7-overlap_bed.f6)
        overlap_bed = overlap_bed.annotate(sv_prop=hl.int(overlap_bed[overlap_field]) / overlap_bed.sv_len, 
                            ref_prop=hl.int(overlap_bed[overlap_field]) / overlap_bed.ref_len)

        # use ref_bed_with_header for annotation column/field names
        ref_bed_with_header = hl.import_table(ref_bed_with_header_uri)                                    
        ref_bed_with_header_idx = range(5, len(fields)-1)
        ref_bed_with_header_mapping = {f"f{ref_bed_with_header_idx[i]}": list(ref_bed_with_header.row)[i].lower().replace(' ', '_') 
                                        for i in range(len(ref_bed_with_header_idx))} | {'sv_prop': f"{annot_name}_overlap"}
        overlap_bed = overlap_bed.rename(ref_bed_with_header_mapping)

        overlap_bed = overlap_bed.annotate(locus=hl.locus(overlap_bed.f0, overlap_bed.f1, genome_build))
        overlap_bed = overlap_bed.annotate(alleles=mt.key_rows_by('locus').rows()[overlap_bed.locus].alleles)
        overlap_bed = overlap_bed.key_by('locus','alleles')

        # annotate original VCF
        annot_fields = list(ref_bed_with_header_mapping.values())[3:]
        mt = mt.annotate_rows(info=mt.info.annotate(
            **{field: overlap_bed[mt.row_key][field] for field in annot_fields}))

        for field in annot_fields:
            header['info'][field] = {'Description': '', 'Number': '.', 'Type': 'String'}
        
        return mt, header

    for (intersect_bed, annot_name, ref_bed_with_header_uri) in zip(intersect_bed_files, annot_names, ref_bed_files):
        mt, header = annotate_with_bed(intersect_bed, annot_name, ref_bed_with_header_uri, mt, header)

    # export annotated VCF
    hl.export_vcf(mt, cohort_prefix + 'annotated.vcf.bgz', metadata=header, tabix=True)
    EOF

    python3 annotate_vcf.py ~{vcf_file} ~{sep=',' intersect_bed_files} ~{bed_files_with_header} \
    ~{cohort_prefix} ~{genome_build} ~{cpu_cores} ~{memory}
    >>>

    output {
        File annotated_vcf = cohort_prefix + 'annotated.vcf.bgz'
        File annotated_vcf_idx = cohort_prefix + 'annotated.vcf.bgz.tbi'
    }
}

task filterVCF {
    input {
        File vcf_file
        File ped_uri
        String genome_build
        String hail_docker
        Float gnomad_af_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    set -eou pipefail
    cat <<EOF > filter_vcf.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    ped_uri = sys.argv[2]
    genome_build = sys.argv[3]
    gnomad_af_threshold = float(sys.argv[4])
    cores = sys.argv[5]
    mem = int(np.floor(float(sys.argv[6])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    def get_transmission(phased_tm):
        phased_tm = phased_tm.annotate_entries(transmission=hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
                hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                            hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                    hl.or_missing(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
        )
        return phased_tm

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    header = hl.get_vcf_metadata(vcf_file)

    # Phasing
    tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
    cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
    tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)
    pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

    tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

    # grab Pathogenic only
    path_tm = phased_tm.filter_rows((phased_tm.info.clinical_interpretation[0].matches('athogenic')) |  # ClinVar P/LP
                (hl.is_defined(phased_tm.info.dbvar_pathogenic[0])) |  # dbVar Pathogenic
                (hl.is_defined(phased_tm.info.gd_sv_name[0])))  # GD region  
    
    path_tm = path_tm.filter_entries((path_tm.proband_entry.GT.is_non_ref()) | 
                                    (path_tm.mother_entry.GT.is_non_ref()) |
                                    (path_tm.father_entry.GT.is_non_ref()))
    path_tm = path_tm.annotate_rows(variant_type='P/LP')
    path_tm = get_transmission(path_tm)
    # Mendel errors
    all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
    path_tm = path_tm.annotate_entries(mendel_code=all_errors_mt[path_tm.row_key, path_tm.col_key].mendel_code)

    # filter
    gnomad_fields = [x for x in list(mt.info) if 'gnomad' in x and 'AF' in x]
    mt = mt.annotate_rows(gnomad_popmax_af=hl.max([mt.info[field] for field in gnomad_fields]))

    filt_mt = mt.filter_rows(((~mt.info.clinical_interpretation[0].matches('enign')) |  # not ClinVar benign
                    (hl.is_missing(mt.info.clinical_interpretation[0]))) &
                (hl.is_missing(mt.info.gnomad_sv_name[0])) &  # not gnomAD benign
                (mt.gnomad_popmax_af <= 0.05))  
    
    # export P/LP TSV
    path_tm.entries().flatten().export(os.path.basename(vcf_file).split('.vcf')[0] + '_path_variants.tsv.gz', delimiter='\t')

    # export filtered VCF
    hl.export_vcf(filt_mt, os.path.basename(vcf_file).split('.vcf')[0] + '.filtered.vcf.bgz', metadata=header, tabix=True)
    EOF
    python3 filter_vcf.py ~{vcf_file} ~{ped_uri} ~{genome_build} ~{gnomad_af_threshold} ~{cpu_cores} ~{memory}
    >>>

    output {
        File pathogenic_tsv = basename(vcf_file, '.vcf.bgz') + '_path_variants.tsv.gz'
        File filtered_vcf = basename(vcf_file, '.vcf.bgz') + '.filtered.vcf.bgz'
        File filtered_vcf_idx = basename(vcf_file, '.vcf.bgz') + '.filtered.vcf.bgz.tbi'
    }
}