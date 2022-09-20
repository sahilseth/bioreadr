

# example -------
# creating mappability files:
# projects/ss_tnbc/ms51/wex/ss_ct_het/cnv/01a_try_purecn.sh

if(FALSE){
    library(pacman)
    p_load(tidyverse)
    p_load(glue, janitor)
    # detach("package:flowr")
    # detach("package:params")
    remotes::install_local("~/Dropbox/public/github_params")
    remotes::install_local("~/Dropbox/public/github_flowr")
    p_load(params, flowr)
    p_load(futile.logger)
    source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/gatkcnv.R")
    opts_flow = new_opts()
    
    bampath = "/rsrch3/home/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam"
    opts_flow$load_toml("~/Dropbox/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml")
    # remotes::install_local("~/Dropbox/public/github_params")
    purecn_exe = opts_flow$get("purecn_exe")
    normaldb_fl = opts_flow$get("purecn_normaldb_fl")
    mappingbias_fl = opts_flow$get("purecn_mappingbias_fl")
    snpblacklist_fl = opts_flow$get("purecn_snpblacklist_fl")
    cosmic_fl = opts_flow$get("purecn_cosmic_fl");cosmic_fl
    intervals = opts_flow$get("purecn_intervals")
    bed = "/rsrch3/home/iacs/sseth/ref/az_ref_beds/ss_downloads/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets_nochr_purecnv_optimized.bed"
    purecn_opts = opts_flow$get("purecn_opts")
    cpu = opts_flow$get("purecn_cpu")
    
    runpath = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_noref/runs"
    trk_full = read_tsv(glue("{runpath}/trk_full.tsv"))
    
    ind = "185_162"
    setwd(glue("{runpath}/{ind}/tmp"))
    wranglr::mkdir("purecn")
    df = trk_full %>% filter(individual == ind)
    df %>% data.frame()
    
    sampleid = df$gatkcnv_prefix_paired %>% basename()
    bam = file.path(bampath, df$bam)
    # tumor_cov = 
    counts_hf5 = gsub("tsv", "hdf5", df$gatkcnv_cnts)
    m2_vcf = glue("{df$mutect2_outprefix_paired}_f2.vcf.gz")
    vcf_fl = "ssm/m1_m2/pcgr/IPCT-S4037-MOON0051-Cap2620-4-NT05___pon_combvcf.pcgr_acmg.grch37.vcf.gz"
    seg_fl = df$gatkcnv_final_seg
    logratio_fl = df$gatkcnv_dn_cr
    
    out_prefix = glue("purecn/{sampleid}")
    purecn_vcf = glue("purecn/{sampleid}_purecn.vcf.gz")
    purecn_cov = glue("purecn/{sampleid}_purecn.cov.loess.txt")
    
}



# single sample purecnv using mutect2 and gatk4

if(FALSE){
    # Internally,createNormalDatabasedetermines the sex of the samplesand trains a PCA that is later used for denoising tumor coverage using Tangent normalization[9]:
    # tumor.coverage.file <- system.file("extdata", "example_tumor.txt", package = "PureCN")
    # tmp = read_tsv(tumor.coverage.file)
    cov = calculateBamCoverageByInterval(bam.file = bam, 
                                         interval.file = intervals, 
                                         output.file = purecn_cov)
    
    # run using gakt counts (guess it might be faster)
    source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/gatkcnv.R")
    cmd <- .collect_read_counts(bam, counts = counts_hf5, gatkcnv_intervals = bed)
    cmd$cmds
    system(glue("module load singularity/3.5.2;{cmd$cmds}"))
    # gc norm:
    debug(PureCN:::.readCoverageGatk4)
    read.length <- 75
    
    G <- h5read("myhdf5file.h5", "df")
    
    x <- read_counts.gatk(tumor_cov) %>% data.frame()
    head(x)
    # intervals <- data.frame(x$intervals$transposed_index_start_end)
    # intervals[, 1] <- x$intervals$indexed_contig_names[intervals[,1] + 1]
    targetCoverage <- GRanges(x[, 1], IRanges(x[, 2], x[, 3]))
    targetCoverage$counts <- x$count
    targetCoverage$coverage <- targetCoverage$counts * read.length * 2
    targetCoverage <- PureCN:::.addAverageCoverage(targetCoverage)
    targetCoverage$on.target <- TRUE
    targetCoverage
    PureCN:::.checkLowCoverage(targetCoverage)
    PureCN:::.checkIntervals(targetCoverage)
    
    # counts = readCoverageFile(counts_hf5)
    out <- correctCoverageBias(targetCoverage, intervals, output.file = purecn_cov, plot.bias = TRUE)
    
    # tumor_cov = read_counts.gatk(tumor_cov)
    # length(normalDB)
    # pool <- calculateTangentNormal(tumor_cov, normalDB)
    
}

purecn.single.m2.gtk4 <- function() {
    
    library(PureCN)
    flog.info("read normal DB")
    
    normalDB <- read_rds(normaldb_fl)
    
    # get coverage -----
    cmd = glue(
        "Rscript {purecn_exe}/Coverage.R ",
        "--outdir purecn ",
        "--bam {bam} ",
        "--intervals {intervals}")
    cmd
    system(cmd)
    
    
    
    # read/re-write segfile
    p_load(dplyr)
    df_seg = read_counts.gatk(seg_fl) %>% 
        mutate(ID = sampleid, 
               contig = switch_chr_names.to_int(contig),
               contig = as.integer(contig)) %>% 
        dplyr::select(ID, chrom = contig, 
                      loc.start = start, loc.end = end, 
                      num.mark = num_points_copy_ratio, 
                      seg.mean = log2_copy_ratio_posterior_50)
    df_seg$chrom
    tail(df_seg)
    # need to write it out, since runabs expects a file :(
    tmpfl = tempfile(fileext = ".txt");tmpfl
    write_tsv(df_seg, tmpfl)
    
    
    # it needs its own cov file
    
    # detect sex
    vcf = readVcf(m2_vcf)
    # rowRanges(vcf)
    # PureCN::getSexFromVcf(vcf)
    # cov2 = "purecn/IPCT-S4037-MOON0051-Cap2620-4-NT05___pon_coverage_loess.txt.gz"
    # unlink(cov2)
    # file.symlink(normalizePath(cov), cov2)
    # PureCN::getSexFromCoverage(cov)
    
    
    # vcf = "IPCT-S4037-MOON0051-Cap2620-4-NT05___pon_combvcf.pcgr_acmg.grch37.pass.vcf.gz"
    vcf = readVcf(vcf_fl)
    df1 = rowRanges(vcf)
    df2 = info(vcf) %>% as_tibble()
    df1
    df1 = head(df2$POPAF) %>% unlist()
    table(df2$DB)
    
    # it does use -log10(popaf)
    # ret$input$centromeres = centromeres
    # ret$input$centromeres
    
    # https://github.com/lima1/PureCN/blob/443408332ca002f7a22b74017319ce2e060417e5/R/filterVcf.R#L460
    p_load(BiocParallel)
    centromeres <- PureCN:::.getCentromerePositions(NULL, "hg19", style = "NCBI")
    BPPARAM <- MulticoreParam(workers = cpu)
    debug(runAbsoluteCN)
    ret <- runAbsoluteCN(
        sampleid = sampleid,
        seg.file = NULL, 
        tumor.coverage.file = tumor,
        interval.file = intervals, 
        vcf.file = vcf_fl, 
        genome = "hg19", 
        
        log.file = "some.log",
        verbose = TRUE, 
        
        max.candidate.solutions = 20, 
        test.purity = seq(0.05, 0.95, by = 0.01),
        post.optimize = TRUE,
        
        centromeres = centromeres,
        
        # germline:read_tsv
        POPAF.info.field = "POPAF",
        min.pop.af = 0.001,
        
        # cosmis details
        cosmic.vcf.file = cosmic_fl,
        args.setPriorVcf = list(min.cosmic.cnt = 6),
        
        
        normalDB = normalDB, 
        plot.cnv = FALSE,
        #    vcf.field.prefix = "PureCN.",
        
        # used for local optima step
        BPPARAM = BPPARAM)
    
    write_rds(ret, file = glue("{out_prefix}.rds"))
    createCurationFile(glue("{out_prefix}.rds"))
    df_summ = purecn_get_summary(ret)
    write_csv(df_summ, glue("{out_prefix}_full.csv"))
    
}


purecn_summary_plots <- function(ret, out_prefix, normaldbfile){
    flog.info("Generating output files...")
    pdf(glue("{out_prefix}.pdf"), width = 10, height = 11)
    plotAbs(ret, type = "all")
    invisible(dev.off())
    
    pdf(glue("{out_prefix}_local_optima.pdf"), width = 5, height = 5)
    plotAbs(ret, 1, type = "overview")
    invisible(dev.off())
    
    seg <- ret$results[[1]]$seg
    seg <- seg[, c(1:6, match("C", colnames(seg)))]
    write.table(seg, file = glue("{out_prefix}_dnacopy.seg"), sep = "\t", quote = FALSE,
                row.names = FALSE)
    
    if (is(ret$results[[1]]$gene.calls, "data.frame")) {
        allAlterations <- callAlterations(ret, all.genes = TRUE)
        
        write.csv(cbind(Sampleid = sampleid, gene.symbol = rownames(allAlterations),
                        allAlterations), row.names = FALSE, file = glue("{out_prefix}_genes.csv"), quote = FALSE)
        if (is.null(normalDB)) 
            normalDB <- readRDS(normaldbfile)
        if (normalDB$version >= 8) {
            file.amps <- glue("{out_prefix}_amplification_pvalues.csv")
            allAmplifications <- callAmplificationsInLowPurity(ret, normalDB,
                                                               all.genes = TRUE, BPPARAM = BPPARAM)
            
            write.csv(cbind(Sampleid = sampleid, gene.symbol = rownames(allAmplifications),
                            allAmplifications), row.names = FALSE, file = file.amps, quote = FALSE)
        }
    } else {
        flog.warn("--intervals does not contain gene symbols. Not generating gene-level calls.")
    }
    
    if (!is.null(ret$input$vcf)) {
        file.vcf <- glue("{out_prefix}_.vcf")
        vcfanno <- predictSomatic(ret, return.vcf = TRUE)
        writeVcf(vcfanno, file = file.vcf)
        bgzip(file.vcf, paste0(file.vcf, ".gz"), overwrite = TRUE)
        indexTabix(paste0(file.vcf, ".gz"), format = "vcf")
        file.csv <- glue("{out_prefix}__variants.csv")
        write.csv(cbind(Sampleid = sampleid, predictSomatic(ret)), file = file.csv,
                  row.names = FALSE, quote = FALSE)
        
        file.loh <- glue("{out_prefix}_loh.csv")
        write.csv(cbind(Sampleid = sampleid, PureCN::callLOH(ret)), file = file.loh,
                  row.names = FALSE, quote = FALSE)
        
        file.pdf <- glue("{out_prefix}_chromosomes.pdf")
        pdf(file.pdf, width = 9, height = 10)
        vcf <- ret$input$vcf[ret$results[[1]]$SNV.posterior$vcf.ids]
        chromosomes <- seqlevelsInUse(vcf)
        chromosomes <- chromosomes[orderSeqlevels(chromosomes)]
        for (chrom in chromosomes) {
            plotAbs(ret, 1, type = "BAF", chr = chrom)
        }
        invisible(dev.off())
    }
    
    # # somatic vcf
    # df_mut = predictSomatic(ret, 2) %>% as_tibble()
    # df_mut %>% write_rds(glue("{out_prefix}_mut.rds"))
    # df_mut %>% filter(gene.symbol == "TP53") %>% data.frame()
    # vcf <- predictSomatic(ret, return.vcf = TRUE)
    # writeVcf(vcf, glue("{out_prefix}_mut.vcf.gz"))
    
    # other metrics
    callableBed <- import(system.file("extdata", "example_callable.bed.gz",package = "PureCN"))
    tmb = callMutationBurden(ret, callable = callableBed)
    tmb2 = callMutationBurden(ret)
    write_tsv(cin, glue("{out_prefix}_cin.tsv"))
    
    # All 4 possible ways to calculate fraction of genome altered
    loh = callLOH(ret, 1) %>% as_tibble()
    write_rds(loh, glue("{out_prefix}_loh.rds"))
    loh
    
    # undebug(PureCN:::.getArmLocations)
    cin <- data.frame(
        cin = callCIN(ret, allele.specific = FALSE, reference.state = "normal"),
        cin.allele.specific = callCIN(ret, reference.state = "normal"),
        cin.ploidy.robust = callCIN(ret, allele.specific = FALSE),
        cin.allele.specific.ploidy.robust = callCIN(ret))
    cin
    write_tsv(cin, glue("{out_prefix}_cin.tsv"))
    
    ret
}

# IPCT-S4037-MOON0051-Cap2619-4-NT28
# 185_136
# 185_162
purecn <- function(
    
    sampleid,
    tumor_cov,
    m2_vcf,
    seg_fl,
    logratio_fl,
    genome,
    
    out_prefix = out_prefix,
    purecn_vcf,
    
    purecn_exe = opts_flow$get("purecn_exe"),
    normaldb_fl = opts_flow$get("purecn_normaldb_fl"),
    mappingbias_fl = opts_flow$get("purecn_mappingbias_fl"),
    snpblacklist_fl = opts_flow$get("purecn_snpblacklist_fl"),
    intervals = opts_flow$get("purecn_intervals"),
    
    purecn_opts = opts_flow$get("purecn_opts"),
    cpu = opts_flow$get("purecn_cpu")
    
    
){
    
    cmd = glue("Rscript {purecn_exe}/PureCN.R ",
               "--out {out_prefix} ",
               "--sampleid {sampleid} ",
               # "--tumor{tumor_cov} ",
               
               "--logratiofile {logratio_fl} ",
               "--segfile {seg_fl} ",
               
               "--mappingbiasfile {mappingbias_fl} ",
               "--vcf {m2_vcf} ",
               "--outvcf {purecn_vcf} ",
               
               "--intervals {intervals} ",
               "--snpblacklist {snpblacklist_fl} ",
               "--genome hg19 ",
               "--funsegmentation Hclust ",
               "--parallel --cores {cpu} ",
               "{purecn_opts}")
    cmd
    system(cmd)
    
}

if(FALSE){
    
    
    
}



