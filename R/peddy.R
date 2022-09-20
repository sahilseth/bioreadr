# https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format

# For these tools, the PED files must contain only the first 6 (mandatory) columns from the PLINK format PED file, and no alleles, like a FAM file in PLINK:
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype

# cat /rsrch3/home/iacs/sseth/flows/SS/tnbc/ms51_wex/qc/peddy/temp.ped
# # familyid indind paternalid maternalid sex phenotype
# FAM001  1  0 0  0  -9


# https://github.com/brentp/peddy
# conda activate peddy

peddy_create_ped <- function(trk, pedfl = glue("{peddy_dir}/wex.ped") %>% tools::file_path_as_absolute()){
        ped = data.frame(familyid = trk$individual,
                    indid = trk$sequencing_id,
                    paternalid = 0,
                    maternalid = 0,
                    sex = 0,
                    phenotype = -9)
    # pedfl = glue("{peddy_dir}/wex.ped") %>% tools::file_path_as_absolute()
    write_tsv(ped, pedfl, col_names = FALSE)
    pedfl
}

peddy <- function(trk, 
    # module load conda_/;source $conda_setup;
    pedfl,
    peddy_exe = "conda activate peddy;peddy",
    peddy_dir = "peddy"){

    # ped="/rsrch3/home/iacs/sseth/flows/SS/tnbc/ms51_wex/qc/peddy/temp.ped"
    # vcf="../gatk_hc/185_355/IPCT-S4011-MOON0051-Cap2492-4-NT58.hc.vcf.gz"

    check_args()
    trk %<>% 
        mutate(peddy_outprefix = file.path(peddy_dir, outprefix), 
                peddy_vcf = paste0(gatk_hc_outprefix, ".hc.g.vcf.gz"),
                gatk_hc_geno_vcf = paste0(gatk_hc_outprefix, ".hc.geno.vcf.gz")) %>% 
        metadata_for_dnaseq()

    i=1
    tmp = lapply(1:nrow(trk), function(i){
        samplename = trk$individual[i]
        vcf = trk$gatk_hc_geno_vcf[i]
        oprefix = trk$peddy_outprefix[i]
        # peddy --plot -p 4 --prefix mystudy $VCF $PED
        cmd = glue("mkdir {peddy_dir};{peddy_exe} --plot -p 4 --prefix {oprefix} {vcf} {pedfl} --loglevel DEBUG")
        cmds <- list(peddy = cmd)
        # mutect_gather_bams = cmd_mutect_gather_bams

        flowmat = to_flowmat(cmds, samplename = samplename) %>%
            dplyr::mutate(cmd = as.character(cmd))

        ret = list(flowmat = flowmat)
        ret
    })
  lst <- purrr::transpose(tmp)
  lst$flowmat %<>% bind_rows()
  lst$trk <- trk

  return(lst)
}
# no test comment

# use this instead:
# https://github.com/brentp/somalier

peddy_read.trk <- function(trk, peddy_dir){
    trk_peddy = trk  %>% mutate(
        peddy_sex = glue("{peddy_dir}/{samplename}/{outprefix}.sex_check.csv"),
        peddy_het = glue("{peddy_dir}/{samplename}/{outprefix}.het_check.csv"),
        peddy_ped = glue("{peddy_dir}/{samplename}/{outprefix}.ped_check.csv"))
    trk_peddy %<>%  tidylog::filter(file.exists(peddy_sex))

    i=3
    i=1
    # flog.threshold(INFO)
    df_peddy = lapply(1:nrow(trk_peddy), function(i){
        flog.debug(trk_peddy$outprefix[i])
        df_sex <- data.table::fread(trk_peddy$peddy_sex[i], data.table = FALSE) %>% 
            dplyr::rename(het_ratio_sex = het_ratio, het_count_sex = het_count)
        df_het <- data.table::fread(trk_peddy$peddy_het[i], data.table = FALSE)
        df = dplyr::full_join(df_sex, df_het, by = "sample_id")
        df
    }) %>% bind_rows() %>% as_tibble() %>% clean_names()
    df_peddy
    
}

somalier <- function(trk){

    cmd1 = glue("somalier extract -d extracted/ --sites sites.vcf.gz -f /data/human/g1k_v37_decoy.fa joint.vcf.gz")

}
