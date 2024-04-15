version 1.0

import "mergeSplitVCF.wdl" as mergeSplitVCF
import "mergeVCFs.wdl" as mergeVCFs
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRareParentsVariantsHail {
    input {
        Array[File]? vep_vcf_files
        Array[File]? vep_annotated_final_vcf
        File lcr_uri
        File ped_uri
        File? meta_uri
        File? trio_uri
        File info_header
        String python_trio_sample_script
        String filter_rare_parents_python_script
        String hail_docker
        String vep_hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Boolean merge_split_vcf
        Boolean sort_after_merge=false
        Float AF_threshold=0.005
        Int AC_threshold=1
        Float csq_af_threshold=0.00001
        Int gq_het_threshold=99
        Int gq_hom_ref_threshold=30
        Int shards_per_chunk=10
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
    }  

    if (!defined(meta_uri)) {
        call makeTrioSampleFiles {
            input:
                python_trio_sample_script=python_trio_sample_script,
                ped_uri=ped_uri,
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker
        }        
    }
    File meta_uri_ = select_first([meta_uri, makeTrioSampleFiles.meta_uri])
    File trio_uri_ = select_first([trio_uri, makeTrioSampleFiles.trio_uri])

    Array[File] vep_files = select_first([vep_vcf_files, vep_annotated_final_vcf])

    if (merge_split_vcf) {
        call mergeSplitVCF.splitFile as splitVEPFiles {
            input:
                file=write_lines(vep_files),
                shards_per_chunk=shards_per_chunk,
                cohort_prefix=cohort_prefix,
                vep_hail_docker=vep_hail_docker
        }
        scatter (chunk_file in splitVEPFiles.chunks) {        
            call mergeVCFs.mergeVCFs as mergeChunk {
                input:
                    vcf_files=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    sort_after_merge=sort_after_merge,
                    runtime_attr_override=runtime_attr_merge_chunk
            }
            call filterRareParentsVariants {
                input:
                    vcf_file=mergeChunk.merged_vcf_file,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    meta_uri=meta_uri_,
                    trio_uri=trio_uri_,
                    info_header=info_header,
                    filter_rare_parents_python_script=filter_rare_parents_python_script,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=cohort_prefix,
                    AC_threshold=AC_threshold,
                    AF_threshold=AF_threshold,
                    csq_af_threshold=csq_af_threshold,
                    gq_het_threshold=gq_het_threshold,
                    gq_hom_ref_threshold=gq_hom_ref_threshold,
                    runtime_attr_override=runtime_attr_filter_vcf
            }
        }
        call helpers.mergeResultsPython as mergeResults_merged {
            input:
                tsvs=filterRareParentsVariants.ultra_rare_parents_tsv,
                hail_docker=hail_docker,
                input_size=size(filterRareParentsVariants.ultra_rare_parents_tsv, 'GB'),
                merged_filename=cohort_prefix+'_ultra_rare_variants.tsv',
                runtime_attr_override=runtime_attr_merge_results
        }
    }
    
    if (!merge_split_vcf) {
        scatter (vcf_file in vep_files) {
            String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
            call filterRareParentsVariants as filterRareParentsVariants_sharded {
                input:
                    vcf_file=vcf_file,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    meta_uri=meta_uri_,
                    trio_uri=trio_uri_,
                    info_header=info_header,
                    filter_rare_parents_python_script=filter_rare_parents_python_script,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=basename(vcf_file, file_ext),
                    AC_threshold=AC_threshold,
                    AF_threshold=AF_threshold,
                    csq_af_threshold=csq_af_threshold,
                    gq_het_threshold=gq_het_threshold,
                    gq_hom_ref_threshold=gq_hom_ref_threshold,
                    runtime_attr_override=runtime_attr_filter_vcf
                    }
        }
        call helpers.mergeResultsPython as mergeResults_sharded {
            input:
                tsvs=filterRareParentsVariants_sharded.ultra_rare_parents_tsv,
                hail_docker=hail_docker,
                input_size=size(filterRareParentsVariants_sharded.ultra_rare_parents_tsv, 'GB'),
                merged_filename=cohort_prefix+'_ultra_rare_variants.tsv',
                runtime_attr_override=runtime_attr_merge_results
        }
    }

    File merged_ultra_rare_parents_tsv = select_first([mergeResults_merged.merged_tsv, mergeResults_sharded.merged_tsv])

    output {
        # File hail_log = filterRareParentsVariants.hail_log
        File ultra_rare_parents_tsv = merged_ultra_rare_parents_tsv
    }
}

task makeTrioSampleFiles {
    input {
        String python_trio_sample_script
        File ped_uri
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    curl ~{python_trio_sample_script} > python_trio_sample_script.py
    python3 python_trio_sample_script.py ~{ped_uri} ~{cohort_prefix} 
    >>>
    
    output {
        File meta_uri = "~{cohort_prefix}_sample_list.txt"
        File trio_uri = "~{cohort_prefix}_trio_list.txt"
        File ped_uri_no_header = "~{cohort_prefix}_no_header.ped"
    }
}

task saveVCFHeader {
    input {
        File vcf_uri
        File info_header
        String sv_base_mini_docker
        Boolean bad_header
    }

    runtime {
        docker: sv_base_mini_docker
    }

    String header_filename = basename(vcf_uri, '.vcf.gz') + '_header.txt'

    command <<<
    bcftools head ~{vcf_uri} > ~{header_filename}
    if [[ "~{bad_header}" == "true" ]]; then
        bcftools head ~{vcf_uri} | grep -v "INFO=" > no_info_header.txt
        cat no_info_header.txt ~{info_header} | sort > ~{header_filename}
    fi
    >>>

    output {
        File header_file = header_filename
    }
}

task filterRareParentsVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File info_header
        String filter_rare_parents_python_script
        String vep_hail_docker
        String cohort_prefix
        Int AC_threshold
        Float AF_threshold
        Float csq_af_threshold
        Int gq_het_threshold
        Int gq_hom_ref_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
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
        docker: vep_hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }


    command <<<
        curl ~{filter_rare_parents_python_script} > filter_rare_variants.py
        python3.9 filter_rare_variants.py ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{AC_threshold} ~{AF_threshold} ~{csq_af_threshold} \
        ~{gq_het_threshold} ~{gq_hom_ref_threshold} > stdout

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
        File ultra_rare_parents_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}
