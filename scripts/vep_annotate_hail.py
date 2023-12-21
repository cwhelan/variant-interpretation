from pyspark.sql import SparkSession
import hail as hl
import sys

vcf_file = sys.argv[1]
vep_annotated_vcf_name = sys.argv[2]
cores = sys.argv[3]
mem = sys.argv[4]

builder = (
                SparkSession 
                .builder
                .config("spark.executor.cores", cores)
                .config("spark.executor.memory", f"{mem}g")
                .enableHiveSupport()
)
                               
spark = builder.getOrCreate()
hl.init(sc=spark.sparkContext)

mt = hl.import_vcf(vcf_file, force_bgz=True, reference_genome='GRCh38')
mt = hl.vep(mt, config='vep_config.json', csq=True)
hl.export_vcf(mt, vep_annotated_vcf_name)