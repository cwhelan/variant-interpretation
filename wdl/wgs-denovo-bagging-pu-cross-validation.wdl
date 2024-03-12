version 1.0

import  "flagRepetitiveRegions.wdl" as flagRepetitiveRegions

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow crossValidateBaggingPU {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_variants_tsv
        File ultra_rare_polyx_vcf
        File repetitive_regions_bed
        Array[String] outlier_samples
        String var_type
        String bagging_pu_source_script
        String cross_validation_script
        String tsv_to_bed_script
        String cohort_prefix
        String sv_base_mini_docker
        String hail_docker
        Boolean remove_regions=true  # remove repetitive regions and multiallelic
        Boolean return_estimators=false
        Array[String]? numeric
        String model_type
        Float vqslod_cutoff=-10
        Float prop_dn=1
        Map[String, String] base_params
        Map[String, String] params
        RuntimeAttr? runtime_attr_bagging_pu
    }

    if (remove_regions) {
        call flagRepetitiveRegions.flagRepetitiveRegions as flagRepetitiveRegions {
            input:
            tsv=ultra_rare_variants_tsv,
            repetitive_regions_bed=repetitive_regions_bed,
            tsv_to_bed_script=tsv_to_bed_script,
            sv_base_mini_docker=sv_base_mini_docker,
            hail_docker=hail_docker
        }

        call removeRegions {
            input:
            tsv=ultra_rare_variants_tsv,
            regions_bed=flagRepetitiveRegions.output_bed,
            hail_docker=hail_docker
        }
    }

    if (!defined(numeric)) {
        Array[String] numeric_default = ['false']
    }
    Array[String] numeric_ = select_first([numeric, numeric_default])

    call crossValidate {
        input:
            vcf_metrics_tsv_final=vcf_metrics_tsv_final,
            ultra_rare_variants_tsv=select_first([removeRegions.regions_removed_tsv, ultra_rare_variants_tsv]),
            ultra_rare_polyx_vcf=ultra_rare_polyx_vcf,
            outlier_samples=outlier_samples,
            var_type=var_type,
            bagging_pu_source_script=bagging_pu_source_script,
            cross_validation_script=cross_validation_script,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            numeric=numeric_,
            vqslod_cutoff=vqslod_cutoff,
            model_type=model_type,
            prop_dn=prop_dn,
            base_params=base_params,
            params=params,
            runtime_attr_override=runtime_attr_bagging_pu
    }

    output {
        File bagging_pu_cv_results = crossValidate.bagging_pu_cv_results
        File bagging_pu_cv_importances = crossValidate.bagging_pu_cv_importances
        File bagging_pu_cv_model = crossValidate.bagging_pu_cv_model
    }
}

task crossValidate {
    input {
        File vcf_metrics_tsv_final
        File ultra_rare_variants_tsv
        File ultra_rare_polyx_vcf
        Array[String] outlier_samples
        Array[String] numeric
        String var_type
        String bagging_pu_source_script
        String cross_validation_script
        String cohort_prefix
        String hail_docker
        Float vqslod_cutoff
        String model_type
        Float prop_dn
        Map[String, String] base_params
        Map[String, String] params
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_metrics_tsv_final, ultra_rare_variants_tsv], "GB")
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
        curl ~{cross_validation_script} > cross_validate.py
        curl ~{bagging_pu_source_script} > baggingPU.py
        python3 cross_validate.py ~{vcf_metrics_tsv_final} ~{ultra_rare_variants_tsv} ~{ultra_rare_polyx_vcf} \
        ~{cohort_prefix} ~{var_type} ~{sep=',' numeric} ~{sep=',' outlier_samples} ~{vqslod_cutoff} ~{model_type} \
        ~{prop_dn} ~{write_map(base_params)} ~{write_map(params)}
    >>>

    output {
        File bagging_pu_cv_results = "~{cohort_prefix}_baggingPU_~{var_type}_~{model_type}_CV_results.tsv"
        File bagging_pu_cv_importances = "~{cohort_prefix}_~{var_type}_~{model_type}_CV_feature_importances.tsv"
        File bagging_pu_cv_model = "~{cohort_prefix}_~{var_type}_~{model_type}_CV_model.pkl"
    }
}

task removeRegions {
    input {
        File tsv
        File regions_bed
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([tsv, regions_bed], "GB")
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
        set -eou pipefail

        cat <<EOF > remove_regions.py 
        import os
        import sys
        import pandas as pd
        import numpy as np

        tsv = sys.argv[1]
        regions_bed = sys.argv[2]

        df = pd.read_csv(tsv, sep='\t')
        regions = pd.read_csv(regions_bed, sep='\t', header=None)

        df['VarKey'] = df[['CHROM','POS','REF','ALT','SAMPLE']].astype(str).agg(':'.join, axis=1)
        df['rm_region'] = df.VarKey.isin(regions[3])

        df['multiallelic'] = (df.DPC_sample!=df.DP_sample)\
                |(df.DPC_mother!=df.DP_mother)\
                |(df.DPC_father!=df.DP_father)

        df = df[(~df.rm_region)&(~df.multiallelic)]
        df.to_csv(f"{os.path.basename(tsv).split('.tsv')[0]}_no_repetitive_regions_multiallelic.tsv", sep='\t', index=False)
        EOF

        python3 remove_regions.py ~{tsv} ~{regions_bed}
    >>>

    output {
        File regions_removed_tsv = basename(tsv, '.tsv') + '_no_repetitive_regions_multiallelic.tsv'
    }

}