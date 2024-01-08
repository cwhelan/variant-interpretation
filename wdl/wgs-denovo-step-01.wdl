version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow step01 {
    input {
        # file can be a list of vcf files or just one vcf file
        File python_trio_sample_script
        File python_preprocess_script
        File flatten_array_script
        File lcr_uri
        File ped_uri
        File info_header
        Array[Array[File]] vep_annotated_final_vcf
        String hail_docker
        String sv_base_mini_docker
        String bucket_id
        String cohort_prefix
        Int? shards_per_chunk
        Boolean bad_header=false
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_flatten
        RuntimeAttr? runtime_attr_preprocess
    }

    call makeTrioSampleFiles {
        input:
            python_trio_sample_script=python_trio_sample_script,
            ped_uri=ped_uri,
            bucket_id=bucket_id,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker
    }

    call flattenArray {
        input:
            nested_array=vep_annotated_final_vcf,
            flatten_array_script=flatten_array_script,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_flatten
    }
    Array[File] vep_annotated_final_vcf_array = flattenArray.flattened_array

    if (defined(shards_per_chunk)) {
        call splitFile {
            input:
                file=write_lines(vep_annotated_final_vcf_array),
                shards_per_chunk=select_first([shards_per_chunk]),
                cohort_prefix=cohort_prefix,
                hail_docker=hail_docker
        }

        scatter (chunk_file in splitFile.chunks) {             
            call mergeVCFs as mergeChunk {
                input:
                    vcf_contigs=read_lines(chunk_file),
                    sv_base_mini_docker=sv_base_mini_docker,
                    cohort_prefix=basename(chunk_file),
                    runtime_attr_override=runtime_attr_merge_vcfs
            }
            call saveVCFHeader as saveVCFHeaderChunk {
                input:
                    vcf_uri=mergeChunk.merged_vcf_file,
                    info_header=info_header,
                    bad_header=bad_header,
                    sv_base_mini_docker=sv_base_mini_docker
            }
            call preprocessVCF as preprocessVCFChunk {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    vcf_uri=mergeChunk.merged_vcf_file,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    hail_docker=hail_docker,
                    header_file=saveVCFHeaderChunk.header_file,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs as mergeChunks {
            input:
                vcf_contigs=preprocessVCFChunk.preprocessed_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_vcfs  
        }
    }

    if (!defined(shards_per_chunk)) {
        scatter (vcf_uri in vep_annotated_final_vcf_array) {
            call saveVCFHeader {
                input:
                    vcf_uri=vcf_uri,
                    info_header=info_header,
                    bad_header=bad_header,
                    sv_base_mini_docker=sv_base_mini_docker
            }
            call preprocessVCF {
                input:
                    python_preprocess_script=python_preprocess_script,
                    lcr_uri=lcr_uri,
                    ped_uri=ped_uri,
                    vcf_uri=vcf_uri,
                    meta_uri=makeTrioSampleFiles.meta_uri,
                    trio_uri=makeTrioSampleFiles.trio_uri,
                    hail_docker=hail_docker,
                    header_file=saveVCFHeader.header_file,
                    runtime_attr_override=runtime_attr_preprocess
            }
        }
        call mergeVCFs {
            input:
                vcf_contigs=preprocessVCF.preprocessed_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix,
                runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    output {
        File meta_uri = makeTrioSampleFiles.meta_uri
        File trio_uri = makeTrioSampleFiles.trio_uri
        File ped_uri_no_header = makeTrioSampleFiles.ped_uri_no_header
        File merged_preprocessed_vcf_file = select_first([mergeChunks.merged_vcf_file, mergeVCFs.merged_vcf_file])
        File merged_preprocessed_vcf_idx = select_first([mergeChunks.merged_vcf_idx, mergeVCFs.merged_vcf_idx])
    }
}

task makeTrioSampleFiles {
    input {
        File python_trio_sample_script
        File ped_uri
        String bucket_id
        String cohort_prefix
        String hail_docker
    }

    runtime {
        docker: hail_docker
    }

    command <<<
    python3 ~{python_trio_sample_script} ~{ped_uri} ~{cohort_prefix} ~{bucket_id}
    >>>
    
    output {
        String meta_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_sample_list.txt"
        String trio_uri = "~{bucket_id}/resources/metadata/~{cohort_prefix}_trio_list.txt"
        String ped_uri_no_header = bucket_id + "/resources/pedigrees/" + cohort_prefix + "_no_header.ped"
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

task preprocessVCF {
    input {
        File python_preprocess_script
        File lcr_uri
        File ped_uri
        File vcf_uri
        File meta_uri
        File trio_uri
        File header_file
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, "GB")
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    python3 ~{python_preprocess_script} ~{lcr_uri} ~{ped_uri} ~{meta_uri} ~{trio_uri} ~{vcf_uri} ~{header_file}
    >>>

    output {
        File preprocessed_vcf = basename(vcf_uri, '.vcf.gz') + '.preprocessed.vcf.bgz'
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_contigs
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    Float input_size = size(vcf_contigs, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
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

    String merged_vcf_name="~{cohort_prefix}.vep.merged.vcf.gz"

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_contigs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat -n --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        # java -jar /usr/picard/picard.jar GatherVcfs -I ~{sep=' -I ' vcf_contigs} -O ~{merged_vcf_name}
        bcftools index -t ~{merged_vcf_name}
    >>>

    output {
        File merged_vcf_file=merged_vcf_name
        File merged_vcf_idx=merged_vcf_name + ".tbi"
    }
}

task flattenArray {
    input {
        Array[Array[File]] nested_array
        File flatten_array_script
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    #  generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    #  generally assume working memory is ~3 * inputs
    #  CleanVcf5.FindRedundantMultiallelics
    
    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: 10,
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        echo ~{nested_array} | python3 ~{flatten_array_script} > flattened.txt
    >>>

    output {
        Array[File] flattened_array = read_lines("flattened.txt")
    }
}

task splitFile {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(file, "GB")
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        split -l ~{shards_per_chunk} ~{file} -a 4 -d "~{cohort_prefix}.shard."
    >>>

    output {
        Array[File] chunks = glob("~{cohort_prefix}.*")
    }
}
