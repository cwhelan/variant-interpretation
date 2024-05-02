version 1.0

## This workflow removes duplicates in a VCF, required for joint genotyping
## Author: asanchis@broadinstitute.org

# WORKFLOW DEFINITION
workflow removeDuplicatesWorkflow {
  input {
    File input_vcf
    String sample_id
    String docker_path
  }

  call removeDuplicatesVCF {
    input:
        input_vcf = input_vcf,
        sample_id = sample_id,
        docker = docker_path
  }

  output {
    File output_vcf = removeDuplicatesVCF.output_vcf
    File output_vcf_index = removeDuplicatesVCF.output_vcf_index
  }
}

# TASK DEFINITIONS

task removeDuplicatesVCF {
  input {

    File input_vcf
    String sample_id
    String docker

    # Runtime parameters
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }

    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

  command {
    set -e

    bcftools norm -D ~{input_vcf} | \
      bcftools norm -m +any | \ ## Merge multiallelic at this stage, pre-joint-genotyping
      bcftools sort -O z -o ~{sample_id}.readgroupadded.uniq.g.vcf.gz

    tabix -p vcf ~{sample_id}.readgroupadded.uniq.g.vcf.gz

  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{sample_id}.readgroupadded.uniq.g.vcf.gz"
    File output_vcf_index = "~{sample_id}.readgroupadded.uniq.g.vcf.gz.tbi"
  }
}