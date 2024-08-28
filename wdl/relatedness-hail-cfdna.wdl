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
        File gnomad_af_resource
        File gnomad_af_resource_idx
        String cohort_prefix
        String relatedness_qc_script
        String plot_relatedness_script
        String sex_qc_script
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        String genome_build
        Int chunk_size=0
        Boolean sort_after_merge=false
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
    }

    Array[File] all_cfdna_vcfs = select_all(cfdna_vcfs)
    scatter(cfdna_vcf in all_cfdna_vcfs) {
        call RemoveRawMutectCalls {
            input:
                vcf = cfdna_vcf,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    Array[File] all_vcfs = select_all(flatten([RemoveRawMutectCalls.out, maternal_vcfs, confirmation_vcfs]))

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
                output_name=prefix + '.somalier.subset.vcf.gz',
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_subset_vcfs
        }
    }

    call mergeVCFs.mergeVCFSamples as mergeSubsetVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            merged_filename=cohort_prefix,
            runtime_attr_override=runtime_attr_merge_vcfs
    }

    call mergeVCFs.mergeVCFSamples as mergeNonSubsetVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            merged_filename=cohort_prefix,
            runtime_attr_override=runtime_attr_merge_vcfs
    }

    File merged_vcf_file = mergeSubsetVCFs.merged_vcf_file

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

    call AnnotateWithGnomadAFs {
        input:
            vcf=mergeNonSubsetVCFs.merged_vcf_file,
            vcf_idx=mergeNonSubsetVCFs.merged_vcf_file_idx,
            gnomad_af_resource=gnomad_af_resource,
            gnomad_af_resource_idx=gnomad_af_resource_idx,
            sv_base_mini_docker=sv_base_mini_docker
    }

    call checkRelatednessRareAlleles {
        input:
            vcf=AnnotateWithGnomadAFs.out,
            vcf_idx=AnnotateWithGnomadAFs.out_idx,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            bucket_id=bucket_id,
            genome_build=genome_build,
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
        File ra_sharing_tsv = checkRelatednessRareAlleles.ra_sharing_tsv
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

task RemoveRawMutectCalls {
    input {
        File vcf
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

    String prefix = basename(vcf, ".vcf.gz")

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
        # mutect calls are the first sample
        bcftools query -l ~{vcf} | awk 'NR > 1' > samples.txt
        bcftools view ~{vcf} -S samples.txt -o ~{prefix}.fm.vcf.gz

    >>>

    output {
        File out = "~{prefix}.fm.vcf.gz"
    }
}

task AnnotateWithGnomadAFs {
    input {
        File vcf
        File vcf_idx
        File gnomad_af_resource
        File gnomad_af_resource_idx
        String sv_base_mini_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3,
                                      disk_gb: 25,
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }

    RuntimeAttr runtime_override = runtime_default
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    String prefix = basename(vcf, ".vcf.gz")

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
        bcftools annotate -a ~{gnomad_af_resource} -c GAF:=AF_popmax -o vcf_annotated_gaf.vcf.gz ~{vcf}
        tabix vcf_annotated_gaf.vcf.gz

    >>>

    output {
        File out = "vcf_annotated_gaf.vcf.gz"
        File out_idx = "vcf_annotated_gaf.vcf.gz.tbi"
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
        File kinship_tsv = cohort_prefix + "_kinship.tsv.gz"
    }
}

task checkRelatednessRareAlleles {
    input {
        File vcf
        File vcf_idx
        String cohort_prefix
        String hail_docker
        String bucket_id
        String genome_build
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
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
        python3 <<CODE
        import numpy as np
        import pysam
        from itertools import combinations

        def compute_rare_allele_sharing(vcf_filename, maf_threshold=0.01):
            # Open the VCF file
            vcf = pysam.VariantFile(vcf_filename)

            # Get list of samples
            samples = list(vcf.header.samples)
            n_individuals = len(samples)

            # Initialize counters for rare alleles
            rare_allele_counts = np.zeros(n_individuals)
            shared_rare_counts = np.zeros((n_individuals, n_individuals))

            # Iterate through variants in the VCF
            for variant in vcf:
                # Check if GAF is present and below threshold
                if 'GAF' in variant.info:
                    gaf = variant.info['GAF'][0]  # Assuming GAF is a single value
                    if gaf < maf_threshold and gaf > 0:
                        genotypes = [variant.samples[sample]['GT'] == (0, 1) or variant.samples[sample]['GT'] == (1, 1) for sample in variant.samples]

                        for i, has_allele in enumerate(genotypes):
                            if has_allele:
                                rare_allele_counts[i] += 1
                                for j in range(i+1, n_individuals):
                                    if genotypes[j]:
                                        shared_rare_counts[i, j] += 1
                                        shared_rare_counts[j, i] += 1

            # Compute pairwise metrics
            results = []
            for (i, j) in combinations(range(n_individuals), 2):
                shared = shared_rare_counts[i, j]
                total_i = rare_allele_counts[i]
                total_j = rare_allele_counts[j]

                prop_a_shared = shared / total_i if total_i > 0 else 0
                prop_b_shared = shared / total_j if total_j > 0 else 0
                jaccard = shared / (total_i + total_j - shared) if (total_i + total_j - shared) > 0 else 0

                results.append((
                    samples[i], samples[j],
                    prop_a_shared, prop_b_shared, jaccard,
                    int(total_i), int(total_j)  # Convert to int for cleaner output
                ))

            return results

        vcf_filename = "~{vcf}"
        sharing_results = compute_rare_allele_sharing(vcf_filename)

        # Write results to TSV file
        output_filename = "~{cohort_prefix}_rare_allele_sharing_results.tsv"
        with open(output_filename, "w") as f:
            f.write("Sample_A\tSample_B\tProp_A_Shared\tProp_B_Shared\tJaccard_Index\tTotal_Rare_A\tTotal_Rare_B\n")
            for result in sharing_results:
                f.write("\t".join(map(str, result)) + "\n")

        print(f"Results written to {output_filename}")
        CODE

    >>>

    output {
        File ra_sharing_tsv = cohort_prefix + "_rare_allele_sharing_results.tsv"
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
