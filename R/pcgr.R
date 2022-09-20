

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

# installation
# ~/Dropbox/misc/privatemodules/seadragon/pcgr_/0.7.0

# nano examples/pcgr_conf.BRCA.toml
# pcgr.py --no-docker --input_vcf examples/tumor_sample.BRCA.vcf.gz  ~/ref/human/pcgr test_out grch37 examples/pcgr_conf.BRCA.toml tumor_sample.BRCA  --force_overwrite


pcgr_colspec_tiers <- function(){
  cols(
    GENOMIC_CHANGE = col_character(),
    GENOME_VERSION = col_character(),
    VCF_SAMPLE_ID = col_character(),
    VARIANT_CLASS = col_character(),
    SYMBOL = col_character(),
    GENE_NAME = col_character(),
    CCDS = col_character(),
    CANONICAL = col_character(),
    ENTREZ_ID = col_double(),
    UNIPROT_ID = col_character(),
    ENSEMBL_TRANSCRIPT_ID = col_character(),
    ENSEMBL_GENE_ID = col_character(),
    REFSEQ_MRNA = col_character(),
    ONCOSCORE = col_double(),
    ONCOGENE = col_logical(),
    TUMOR_SUPPRESSOR = col_logical(),
    DISGENET_CUI = col_character(),
    DISGENET_TERMS = col_character(),
    CONSEQUENCE = col_character(),
    PROTEIN_CHANGE = col_character(),
    PROTEIN_DOMAIN = col_character(),
    CDS_CHANGE = col_character(),
    HGVSp = col_character(),
    HGVSc = col_character(),
    EFFECT_PREDICTIONS = col_character(),
    MUTATION_HOTSPOT = col_character(),
    MUTATION_HOTSPOT_TRANSCRIPT = col_character(),
    MUTATION_HOTSPOT_CANCERTYPE = col_character(),
    INTOGEN_DRIVER_MUT = col_logical(),
    VEP_ALL_CONSEQUENCE = col_character(),
    DBSNPRSID = col_character(),
    COSMIC_MUTATION_ID = col_character(),
    TCGA_PANCANCER_COUNT = col_double(),
    TCGA_FREQUENCY = col_character(),
    ICGC_PCAWG_OCCURRENCE = col_character(),
    CHEMBL_COMPOUND_ID = col_character(),
    CHEMBL_COMPOUND_TERMS = col_character(),
    CLINVAR = col_character(),
    CLINVAR_CLNSIG = col_character(),
    GLOBAL_AF_GNOMAD = col_double(),
    GLOBAL_AF_1KG = col_double(),
    CALL_CONFIDENCE = col_logical(),
    DP_TUMOR = col_double(),
    AF_TUMOR = col_double(),
    DP_NORMAL = col_double(),
    AF_NORMAL = col_double(),
    TIER = col_character(),
    TIER_DESCRIPTION = col_character(),
    TDP = col_double(),
    NDP = col_double(),
    TAF = col_double(),
    NAF = col_double(),
    TAF_FWD = col_double(),
    TAF_REV = col_double(),
    set = col_character()
  )
}

pcgr <- function(vcf = tmp$vcf_out_combvcf, 
                 outprefix_vcf = tmp$outprefix_vcf,
                 pcgr_vcf = tmp$pcgr_vcf,
                 pcgr_tsv = tmp$pcgr_tsv,
                 pcgr_log = tmp$pcgr_log,
                 pcgr_dir = tmp$pcgr_dir,
                 conf = "~/Dropbox/public/flowr/my.ultraseq/pipelines/dnaseq/ss_cl_het/pcgr_conf.BRCA_matched_merged.toml",
                 pcgr = "/rsrch3/home/iacs/sseth/.conda/envs/pcgr_py36/bin/pcgr.py"){
  
  # usage: pcgr.py [options] <PCGR_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY> <CONFIG_FILE> <SAMPLE_ID>
  # positional arguments:
  # pcgr_dir              PCGR base directory with accompanying data directory,
  # output_dir            Output directory
  # {grch37,grch38}       Genome assembly build: grch37 or grch38
  # configuration_file    PCGR configuration file (TOML format)
  # sample_id             Tumor sample/cancer genome identifier - prefix for output files
  # samplename = "185_145_T0-D"
  # samplename = "185_145_T1-D"
  # vcf = glue("variants/{samplename}.vcf.gz")
  
  
  # seq of files (from PCGR) --> ready vcf (temp) ---> maf
  # pcgr_acmg.grch37.vcf.gz
  # pcgr_acmg.grch37_pass.vcf.gz -> pcgr_acmg.grch37_pass.tsv.gz
  # pcgr_acmg.grch37.snvs_indels.tiers.tsv -> json > html
  flog.info("pcgr: running pcgr")
  # pcgr_mod <- "module load python;source activate pcgr_py36"
  pcgr_mod = "module load conda_/3;source $conda_setup;conda activate pcgr_py36"
  cmd = glue("{pcgr_mod};pcgr.py --no-docker --input_vcf {vcf} $HOME/ref/human/pcgr ",
             "{pcgr_dir} grch37 {conf} {outprefix_vcf} --force_overwrite > {pcgr_log}") #| tee
  flog.debug(cmd)
  system(cmd)
  
  flog.info("PCGR: create TSV")
  # tmp = gsub(".gz$", "", tmp$tsv_pcgr_snp)
  tsv = tools::file_path_sans_ext(pcgr_tsv)
  cmd=glue("{pcgr_mod};vcf2tsv.py --keep_rejected_calls --compress --print_data_type_header ",
           "{pcgr_vcf} {tsv} >> {pcgr_log}")
  flog.debug(cmd)
  system(cmd)
}
