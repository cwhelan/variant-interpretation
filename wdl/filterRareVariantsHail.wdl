version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterRareVariantsHail {
    input {
        Array[Array[File]] vep_annotated_final_vcf
        Array[Array[File]] vep_annotated_final_vcf_idx
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File filter_rare_variants_python_script
        File info_header
        String vep_hail_docker
        String sv_base_mini_docker
        String cohort_prefix
        Boolean merge_annotated_vcfs=false
        Boolean bad_header=false
        Float AF_threshold=0.005
        Int AC_threshold=2
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_merge_results
    }

    Array[Pair[Array[File], Array[File]]] vep_annotated_final = zip(vep_annotated_final_vcf, vep_annotated_final_vcf_idx)
    
    if (merge_annotated_vcfs) {
        scatter (vep_annotated_shard in vep_annotated_final) {
            Array[File] cohort_vcf_files = vep_annotated_shard.left
            Array[File] cohort_vcf_idx = vep_annotated_shard.right
            call mergeVCFs as mergeSharded {
                input:
                    vcf_files=cohort_vcf_files,
                    vcf_files_idx=cohort_vcf_idx,
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=cohort_prefix,
                    merge_or_concat='concat',
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
        }

        if (length(mergeSharded.merged_vcf_file) > 1) {
            call mergeVCFs as mergeCohort {
            input:
                vcf_files=mergeSharded.merged_vcf_file,
                vcf_files_idx=mergeSharded.merged_vcf_idx,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                merge_or_concat='concat',
                runtime_attr_override=runtime_attr_merge_vcfs
            }
        }
        
        call filterRareVariants {
            input:
                vcf_file=select_first([mergeCohort.merged_vcf_file, mergeSharded.merged_vcf_file[0]]),
                lcr_uri=lcr_uri,
                ped_uri=ped_uri,
                meta_uri=meta_uri,
                trio_uri=trio_uri,
                info_header=info_header,
                filter_rare_variants_python_script=filter_rare_variants_python_script,
                vep_hail_docker=vep_hail_docker,
                cohort_prefix=cohort_prefix,
                AC_threshold=AC_threshold,
                AF_threshold=AF_threshold,
                bad_header=bad_header,
                runtime_attr_override=runtime_attr_filter_vcf
        }
    }

    if (!merge_annotated_vcfs) {
        scatter (vep_annotated_shard in vep_annotated_final) {
            Array[File] cohort_vcf_files_ = vep_annotated_shard.left
            Array[File] cohort_vcf_idx_ = vep_annotated_shard.right
            scatter (vcf_file in cohort_vcf_files_) {
                call filterRareVariants as filterRareVariants_sharded {
                    input:
                        vcf_file=vcf_file,
                        lcr_uri=lcr_uri,
                        ped_uri=ped_uri,
                        meta_uri=meta_uri,
                        trio_uri=trio_uri,
                        info_header=info_header,
                        filter_rare_variants_python_script=filter_rare_variants_python_script,
                        vep_hail_docker=vep_hail_docker,
                        cohort_prefix=basename(vcf_file),
                        AC_threshold=AC_threshold,
                        AF_threshold=AF_threshold,
                        bad_header=bad_header,
                        runtime_attr_override=runtime_attr_filter_vcf
                }
            }

            call mergeResults as mergeShardResults {
                input:
                    ultra_rare_variants_tsvs=filterRareVariants_sharded.ultra_rare_variants_tsv,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=cohort_prefix,
                    runtime_attr_override=runtime_attr_merge_results
            }
        }
        if (length(mergeShardResults.ultra_rare_variants_tsv) > 1) {
            call mergeResults {
                input: 
                    ultra_rare_variants_tsvs=mergeShardResults.ultra_rare_variants_tsv,
                    vep_hail_docker=vep_hail_docker,
                    cohort_prefix=cohort_prefix,
                    runtime_attr_override=runtime_attr_merge_results
            }
        }
    }

    output {
        # File hail_log = filterRareVariants.hail_log
        File ultra_rare_variants_tsv = select_first([filterRareVariants.ultra_rare_variants_tsv, mergeResults.ultra_rare_variants_tsv, mergeShardResults.ultra_rare_variants_tsv])
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

task mergeVCFs {
    input {
        Array[File] vcf_files
        Array[File] vcf_files_idx
        String sv_base_mini_docker
        String cohort_prefix
        String merge_or_concat    
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_files, "GB")
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
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"

    String merge_or_concat_new = if merge_or_concat == 'concat' then 'concat -n'  else merge_or_concat

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools ~{merge_or_concat_new} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}

task filterRareVariants {
    input {
        File vcf_file
        File lcr_uri
        File ped_uri
        File meta_uri
        File trio_uri
        File info_header
        File filter_rare_variants_python_script
        String vep_hail_docker
        String cohort_prefix
        Int AC_threshold
        Float AF_threshold
        Boolean bad_header
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

    String header_filename = basename(vcf_file, '.vcf.gz') + '_header.txt'

    command <<<
        /opt/vep/bcftools/bcftools head ~{vcf_file} > ~{header_filename}
        if [[ "~{bad_header}" == "true" ]]; then
            /opt/vep/bcftools/bcftools head ~{vcf_file} | grep -v "INFO=" > no_info_header.txt
            cat no_info_header.txt ~{info_header} | LC_ALL=C sort > ~{header_filename}
        fi

        python3.9 ~{filter_rare_variants_python_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_file} \
        ~{cohort_prefix} ~{cpu_cores} ~{memory} ~{AC_threshold} ~{AF_threshold} ~{header_filename}

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File hail_log = "hail_log.txt"
        File ultra_rare_variants_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}

task mergeResults {
    input {
        Array[File] ultra_rare_variants_tsvs
        String vep_hail_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ultra_rare_variants_tsvs, "GB")
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
        head -n 1 ~{ultra_rare_variants_tsvs[0]} > "~{cohort_prefix}_ultra_rare_variants.tsv"; 
        tail -n +2 -q ~{sep=' ' ultra_rare_variants_tsvs} >> "~{cohort_prefix}_ultra_rare_variants.tsv"
    >>>

    output {
        File ultra_rare_variants_tsv = cohort_prefix + '_ultra_rare_variants.tsv'
    }
}