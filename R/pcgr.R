

# conda load pcgr_py36
# 

# conda activate pvactools

# wget --no-check-certificate https://drive.google.com/uc?export=download&confirm=WFJt&id=1OL5C994HDaeadASz7KzMhPoXfdSiyhNy -O pcgr.databundle.grch37.20181119.tgz
# 403  wget --no-check-certificate 'https://drive.google.com/uc?export=download&confirm=WFJt&id=1OL5C994HDaeadASz7KzMhPoXfdSiyhNy' -O pcgr.databundle.grch37.20181119.tgz
# 407  tar -zxvf pcgr.databundle.grch37.20181119.tgz
# 409  tar -zxvf pcgr.databundle.grch38.20181119.tgz
# 425  gzip -dc pcgr.databundle.grch37.20181119.tgz | tar xvf -


# seadragon test
# module load python
# source activate pcgr_py36

# nano examples/pcgr_conf.BRCA.toml
# pcgr.py --no-docker --input_vcf examples/tumor_sample.BRCA.vcf.gz     ~/ref/human/pcgr test_out grch37 examples/pcgr_conf.BRCA.toml tumor_sample.BRCA  --force_overwrite



pcgr <- function(vcf, 
                 conf = "/rsrch3/home/iacs/sseth/apps/pcgr/0.7.0/examples/pcgr_conf.BRCA.toml"){
  
  # usage: pcgr.py [options] <PCGR_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CONFIG_FILE> <SAMPLE_ID>
  # positional arguments:
  # pcgr_dir              PCGR base directory with accompanying data directory,
  # output_dir            Output directory
  # {grch37,grch38}       Genome assembly build: grch37 or grch38
  # configuration_file    PCGR configuration file (TOML format)
  # sample_id             Tumor sample/cancer genome identifier - prefix for output files
  # samplename = "185_145_T0-D"
  # samplename = "185_145_T1-D"
  vcf = glue("variants/{samplename}.vcf.gz")
  cmd = glue("pcgr.py --no-docker --input_vcf {vcf} ", 
             "$HOME/ref/human/pcgr variants grch37 {conf} {samplename}_pcgr  --force_overwrite")
  cmd
  
  
}