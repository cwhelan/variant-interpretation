version 1.0

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

workflow Relatedness {
    input {
        Array[File?] cfdna_vcfs
        Array[File?] maternal_vcfs
        Array[File?] confirmation_vcfs
        File ped_uri
        File sites_uri
        File bed_file
        String cohort_prefix
        String relatedness_qc_script
        String plot_relatedness_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        String genome_build
        Int chunk_size=0
        Boolean sort_after_merge=false
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
    }

    Array[File] all_vcfs = select_all(flatten([cfdna_vcfs, maternal_vcfs, confirmation_vcfs]))

    call DedupVcfs {
        input:
            vcf_list = all_vcfs,
            sv_base_mini_docker = sv_base_mini_docker
    }

    scatter (vcf_uri in DedupVcfs.vcfs) {
        String filename = basename(vcf_uri)
        String prefix = if (sub(filename, ".gz", "")!=filename) then basename(filename, ".vcf.gz") else basename(filename, ".vcf.bgz")
        call helpers.subsetVCFs as subsetVCFs {
            input:
                bed_file=bed_file,
                vcf_uri=vcf_uri,
                vcf_idx=vcf_uri+'.tbi',
                output_name=prefix + '.somalier.subset.vcf.gz',
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_subset_vcfs
        }
    }

    call mergeVCFs.mergeVCFSamples as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            merged_filename=cohort_prefix,
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_vcfs
    }

    File merged_vcf_file = mergeVCFs.merged_vcf_file

    call checkRelatedness {
        input:
        vcf_uri=merged_vcf_file,
        sites_uri=sites_uri,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        relatedness_qc_script=relatedness_qc_script,
        hail_docker=hail_docker,
        bucket_id=bucket_id,
        genome_build=genome_build,
        score_table=false,
        runtime_attr_override=runtime_attr_check_relatedness
    }

    call plotRelatedness {
        input:
        kinship_tsv=checkRelatedness.kinship_tsv,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker,
        chunk_size=chunk_size,
        runtime_attr_override=runtime_attr_plot_relatedness
    }

    output {
        File relatedness_qc = checkRelatedness.relatedness_qc
        File kinship_tsv = checkRelatedness.kinship_tsv
        File relatedness_plot = plotRelatedness.relatedness_plot
    }
}

task DedupVcfs {
    input {
        Array[String] vcf_list
        String sv_base_mini_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3,
                                      disk_gb: 10,
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }

    RuntimeAttr runtime_override = runtime_default
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        sort < ~{write_lines(vcf_list)} | uniq > vcf.dedup.list
    >>>

    output {
        Array[File] vcfs = read_lines("vcf.dedup.list")
    }
}

task checkRelatedness {
    input {
        File vcf_uri
        File ped_uri
        File sites_uri
        String cohort_prefix
        String relatedness_qc_script
        String hail_docker
        String bucket_id
        String genome_build
        String score_table=false
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_uri, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0

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
        curl ~{relatedness_qc_script} > check_relatedness.py
        python3 check_relatedness.py ~{vcf_uri} ~{sites_uri} ~{cohort_prefix} ~{ped_uri} ~{cpu_cores} ~{memory} \
        ~{bucket_id} ~{score_table} ~{genome_build} > stdout
    >>>

    output {
        File relatedness_qc = cohort_prefix + "_relatedness_qc.ped"
        File kinship_tsv = cohort_prefix + "_kinship.tsv"
    }
}

task plotRelatedness {
    input {
        File kinship_tsv
        File ped_uri
        String cohort_prefix
        String plot_relatedness_script
        String hail_docker
        Int chunk_size
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(ped_uri, "GB")
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
        curl ~{plot_relatedness_script} > plot_relatedness.py
        python3 plot_relatedness.py ~{kinship_tsv} ~{cohort_prefix} ~{ped_uri} ~{chunk_size} > stdout
    >>>

    output {
        File relatedness_plot = "~{cohort_prefix}_relatedness_ibd0_kinship.png"
    }
}
