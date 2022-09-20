bioc_genomeinfo_create_cache <- function(){
    # final solution
    # create new env -----
    library(pacman)
    p_load(GenomeInfoDb)
    ls(GenomeInfoDb:::.UCSC_cached_chrom_info)

    # ncbi
    Seqinfo(genome = "GRCh37")
    Seqinfo(genome = "GRCh38")
    Seqinfo(genome = "GRCh38.p12")
    Seqinfo(genome = "hg18")
    Seqinfo(genome = "hg19")
    Seqinfo(genome = "hg38")
    Seqinfo(genome = "GCF_000001405.25")
    human_chrom_info <- getChromInfoFromEnsembl("hsapiens")


    setwd("/rsrch3/home/iacs/sseth/R/bioc_cache")
    save(
        list = ls(all = TRUE, envir = GenomeInfoDb:::.UCSC_cached_chrom_info),
        envir = GenomeInfoDb:::.UCSC_cached_chrom_info,
        file = "GenomeInfoDb___UCSC_cached_chrom_info.rds"
    )
    save(
        list = ls(all = TRUE, envir = GenomeInfoDb:::.NCBI_cached_chrom_info),
        envir = GenomeInfoDb:::.NCBI_cached_chrom_info,
        file = "GenomeInfoDb___NCBI_cached_chrom_info.rds"
    )
    save(
        list = ls(all = TRUE, envir = GenomeInfoDb:::.ENSEMBL_cached_chrom_info),
        envir = GenomeInfoDb:::.ENSEMBL_cached_chrom_info,
        file = "GenomeInfoDb___ENSEMBL_cached_chrom_info.rds"
    )
    save(
        list = ls(all = TRUE, envir = GenomeInfoDb:::.Ensembl_FTP_cached_core_dbs),
        envir = GenomeInfoDb:::.Ensembl_FTP_cached_core_dbs,
        file = "GenomeInfoDb___Ensembl_FTP_cached_core_dbs.rds"
    )
    save(
        list = ls(all = TRUE, envir = GenomeInfoDb:::.Ensembl_FTP_cached_releases),
        envir = GenomeInfoDb:::.Ensembl_FTP_cached_releases,
        file = "GenomeInfoDb___Ensembl_FTP_cached_releases.rds"
    )

}

bioc_genomeinfo_load_cache <- function(){
    p_load(GenomeInfoDb)
    load("~/R/bioc_cache/GenomeInfoDb___UCSC_cached_chrom_info.rds",
        envir = GenomeInfoDb:::.UCSC_cached_chrom_info, verbose = TRUE
    )
    load("~/R/bioc_cache/GenomeInfoDb___NCBI_cached_chrom_info.rds",
        envir = GenomeInfoDb:::.NCBI_cached_chrom_info, verbose = TRUE
    )
    load("~/R/bioc_cache/GenomeInfoDb___ENSEMBL_cached_chrom_info.rds",
        envir = GenomeInfoDb:::.ENSEMBL_cached_chrom_info, verbose = TRUE
    )
    load("~/R/bioc_cache/GenomeInfoDb___Ensembl_FTP_cached_core_dbs.rds",
        envir = GenomeInfoDb:::.Ensembl_FTP_cached_core_dbs, verbose = TRUE
    )
    # ls(GenomeInfoDb:::.UCSC_cached_chrom_info)
    # ls(GenomeInfoDb:::.NCBI_cached_chrom_info)
    Seqinfo(genome = "GRCh37")

}