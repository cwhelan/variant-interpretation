version 1.0 

workflow step4 {
    input {
        File ped_uri
        Array[Array[File]] split_trio_vcfs
        String trio_denovo_docker
    }
    scatter (vcf_files in split_trio_vcfs) {
        scatter (vcf_file in vcf_files) {
            call trio_denovo {
                input:
                    ped_uri=ped_uri,
                    vcf_file=vcf_file,
                    trio_denovo_docker=trio_denovo_docker
            }
        }
        call combineOutputVCFs {
            input:
                out_vcfs=trio_denovo.out_vcf
        }
    }

    output {
        Array[Array[File]] trio_denovo_vcf = combineOutputVCFs.trio_denovo_vcf
    }
}

task trio_denovo {
    input {
        File ped_uri
        File vcf_file
        String trio_denovo_docker
    }

    runtime {
        docker: trio_denovo_docker
    }

    command {
        /src/wgs_denovo/triodenovo/triodenovo-fix/src/triodenovo --ped ~{ped_uri} --in_vcf ~{vcf_file} --out_vcf ~{basename(vcf_file, '.vcf') + '_trio_denovo.vcf'}
        bgzip ~{basename(vcf_file, '.vcf') + '_trio_denovo.vcf'}
    }

    output {
        File out_vcf = basename(vcf_file, '.vcf') + '_trio_denovo.vcf.gz'
    }
}

task combineOutputVCFs {
    input {
        Array[File] out_vcfs
    }

    command {
        mkdir -p tmp_out_vcfs
        mv ~{out_vcfs} tmp_out_vcfs/
    }

    output {
        Array[File] trio_denovo_vcf = glob('tmp_out_vcfs/*')
    }
}
