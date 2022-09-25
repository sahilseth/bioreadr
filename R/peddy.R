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


