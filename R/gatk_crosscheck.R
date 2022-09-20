

# create vcfs for all bams 
gatk_fingerprint <- function(trk,
                            gatk_fp_dir = "gatk_fingerprint",
                            ref_fasta = opts_flow$get('ref_fasta'),
                            java_exe = "/risapps/noarch/jdk/jdk1.8.0_131/bin/java",
                            picard_jar = "/rsrch3/home/iacs/sseth/apps/picard/2.23.4/picard.jar",
                            hapmap_db = "~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map"

){

    trk %<>% 
        mutate(gatk_fp_outprefix = file.path(gatk_fp_dir, outprefix), 
                gatk_fp_vcf = paste0(gatk_fp_outprefix, ".fp.vcf.gz")) %>% 
        metadata_for_dnaseq()

    tmp = lapply(1:nrow(trk), function(i){
    samplename = trk$samplename[i]
    bam = trk$bam[i]
    vcf = trk$gatk_fp_vcf[i]
    cmd = glue("mkdir {gatk_fp_dir};{java_exe} -jar {picard_jar}  ExtractFingerprint -I {bam} -H {hapmap_db} -R {ref_fasta} -O {vcf}")
    # system(cmd)
    flowmat = to_flowmat(list(gatk_fp = cmd), samplename = samplename) %>%
        dplyr::mutate(cmd = as.character(cmd))
    ret = list(flowmat = flowmat)
    })
    lst = purrr::transpose(tmp)
    lst$flowmat %<>% bind_rows()
    lst$trk = trk

    return(lst)
}

# https://gatk.broadinstitute.org/hc/en-us/articles/360040098592-CrosscheckFingerprints-Picard-
# This tool calculates the LOD score for identity check between "groups" of data 
# in the INPUT files as defined by the CROSSCHECK_BY argument. 
# A positive value indicates that the data seems to have come from the same 
# individual or, in other words the identity checks out. 
# The scale is logarithmic (base 10), so a LOD of 6 indicates that it is 1,000,000 
# more likely that the data matches the genotypes than not. 
# A negative value indicates that the data do not match. 
# A score that is near zero is inconclusive and can result from low coverage 
# or non-informative genotypes. Each group is assigned a sample identifier 
# (for SAM this is taken from the SM tag in the appropriate readgroup header line, 
# for VCF this is taken from the column label in the file-header. 
# After combining all the data from the same group together, 
# an all-against-all comparison is performed. 
# Results are categorized as one of EXPECTED_MATCH, EXPECTED_MISMATCH, 
# UNEXPECTED_MATCH, UNEXPECTED_MISMATCH, or AMBIGUOUS depending on the 
# LOD score and on whether the sample identifiers of the groups agree: 
# LOD scores that are less than LOD_THRESHOLD are considered mismatches, 
# and those greater than -LOD_THRESHOLD are matches (between is ambiguous). 
# If the sample identifiers are equal, the groups are expected to match. 
# They are expected to mismatch otherwise. 
# The identity check makes use of haplotype blocks defined in the 
# HAPLOTYPE_MAP file to enable it to have higher statistical power for 
# detecting identity or swap by aggregating data from several SNPs in the 
# haplotype block. This enables an identity check of samples with very 
# low coverage (e.g. ~1x mean coverage). ascat
# When provided a VCF, the identity check looks at the PL, GL and GT fields 
# (in that order) and uses the first one that it finds.
gatk_crosscheck <- function(fls, 
                            outdir,
                            gatk4_exe = opts_flow$get("gatk4_exe"),
                            gatk4_sif = opts_flow$get("gatk4_sif"),
                            hapmap_db = "~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map", 
                            cores = 1){

    check_args()
    gatk_crosscheck_opts = glue("--LOD_THRESHOLD=-5 --EXPECT_ALL_GROUPS_TO_MATCH=true --CROSSCHECK_BY SAMPLE --NUM_THREADS {cores}")
    # singularity exec --bind /risapps --bind /rsrch3 /rsrch3/home/iacs/sseth/sifs/gatk_4.1.7.0.sif gatk CrosscheckFingerprints \
    # --INPUT=tmp/IPCT-S2006-MOON0051-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG.bwa_recalibed.bam \
    # --INPUT=tmp/IPCT-S4008-MOON0051-Cap2043-4-ID01_190415-A00422-0045-AHK77HDSXX-4-ATCACG.bwa_recalibed.bam \
    # --HAPLOTYPE_MAP=~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map --LOD_THRESHOLD=-5  --EXPECT_ALL_GROUPS_TO_MATCH=true --OUTPUT=tmp/sample.crosscheck_metrics
    # bams = trk$bam
    fl_list = paste("-I ", fls, sep = "", collapse = " ")
    hapmap_db = tools::file_path_as_absolute(hapmap_db)
    # with lots of files, we may hit the singularity limit
    cmd1 = glue("module load singularity/3.5.2;{gatk4_exe} --java-options '-Xmx32g' CrosscheckFingerprints {fl_list} --OUTPUT {outdir}/gatk_crosscheck.metrics --HAPLOTYPE_MAP={hapmap_db} {gatk_crosscheck_opts}")
    # incase files are too many
    cmd2 = glue("gatk --java-options '-Xmx32g' CrosscheckFingerprints {fl_list} --OUTPUT {outdir}/gatk_crosscheck.metrics --HAPLOTYPE_MAP={hapmap_db} {gatk_crosscheck_opts}")
        # echo "long command (>120k chars)" | singularity exec image bash -c "read -u 0 line; set -- \$line; exec \"\$@\""
    cmd2_wrap = glue("module load singularity/3.5.2;singularity exec --bind /risapps --bind /rsrch3 {gatk4_sif} {outdir}/gatk_crosscheck.sh ")

    # system(cmd)
    list(cmd1 = cmd1, cmd2 = cmd2, cmd2_wrap = cmd2_wrap)
}

# not working at the moment
gatk_crosscheck2 <- function(bams,
                            odir = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/fingerprints",
                            gatk4_exe = opts_flow$get("gatk4_exe"),
                            hapmap_db = "/rsrch3/home/iacs/sseth/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.interval_list",
                            ref_fasta = opts_flow$envir$ref_fasta,
                            cores = 1) {
    check_args()
    # first create vcfs:
    # Single-sample GVCF calling with allele-specific annotations
    # bams = head(bams)
    fingerprint_vcfs = basename(bams) %>% tools::file_path_sans_ext() %>%
        file.path(odir, .) %>% paste0(".vcf")
    fingerprint_vcfs

    # https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035890551-Allele-specific-annotation-and-filtering-of-germline-short-variants
    hapmap_db = "~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map"
    cmd_splt = glue("module load singularity/3.5.2;{gatk4_exe} --java-options '-Xmx32g' HaplotypeCaller -L {hapmap_db} -I {bams} -O {fingerprint_vcfs} ",
        "-R {ref_fasta} ",
         "-ERC GVCF --native-pair-hmm-threads 4 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation") # -G Standard -G AS_Standard 
    system(cmd_splt[1])

    gatk_crosscheck_opts <- glue("--LOD_THRESHOLD=-5 --EXPECT_ALL_GROUPS_TO_MATCH=true --CROSSCHECK_BY SAMPLE --NUM_THREADS {cores}")
    # singularity exec --bind /risapps --bind /rsrch3 /rsrch3/home/iacs/sseth/sifs/gatk_4.1.7.0.sif gatk CrosscheckFingerprints \
    # --INPUT=tmp/IPCT-S2006-MOON0051-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG.bwa_recalibed.bam \
    # --INPUT=tmp/IPCT-S4008-MOON0051-Cap2043-4-ID01_190415-A00422-0045-AHK77HDSXX-4-ATCACG.bwa_recalibed.bam \
    # --HAPLOTYPE_MAP=~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map --LOD_THRESHOLD=-5  --EXPECT_ALL_GROUPS_TO_MATCH=true --OUTPUT=tmp/sample.crosscheck_metrics
    # bams = trk$bam
    bam_list <- paste("--INPUT ", bams, sep = "", collapse = " ")
    hapmap_db <- tools::file_path_as_absolute(hapmap_db)
    cmd <- glue("module load singularity/3.5.2;{gatk4_exe} --java-options '-Xmx32g' CrosscheckFingerprints {bam_list} --OUTPUT gatk_crosscheck.metrics --HAPLOTYPE_MAP={hapmap_db} {gatk_crosscheck_opts}")
    # system(cmd)
}


gatk_crosscheck_ann <- function(df_cross, df_ann){
    source("~/Dropbox/public/github_wranglr/R/add_prefix_suffix_colnames.R")
    wranglr::expect_columns(df_ann, c("individual", "sample_type", "sequencing_id"))

    df_cross_ann <- add_prefix_suffix_column(df_ann, prefix = "left_") %>%
        tidylog::left_join(df_cross, ., by = c("left_group_value" = "left_sequencing_id"))
    df_cross_ann <- add_prefix_suffix_column(df_ann, prefix = "right_") %>%
        tidylog::left_join(df_cross_ann, ., by = c("right_group_value" = "right_sequencing_id"))
}

gatk_crosscheck_summary <- function(df_cross_ann, outdir){
    # unexpected mismatches --------
    df_cross_nomatch = df_cross_ann %>%
        tidylog::filter(left_individual == right_individual) %>%
        tidylog::filter(lod_score < 0) %>% 
        arrange(left_individual, left_sample_type, 
                right_individual, right_sample_type)

    df_cross_mis = df_cross_ann %>%
        tidylog::filter(lod_score > 0) %>%
        tidylog::filter(left_individual != right_individual) %>% 
        arrange(left_individual, left_sample_type, 
                right_individual, right_sample_type)

    df_cross_match = df_cross_ann %>%
        tidylog::filter(lod_score > 0 | left_individual == right_individual)

    source("~/Dropbox/public/github_params/R/write_sheet.R")
    p_load(openxlsx)
    ret = list(df_cross_nomatch = df_cross_nomatch, df_cross_mis = df_cross_mis, df_cross_match = df_cross_match)
    write_sheets(ret, glue("{outdir}/gatk_crosscheck.xlsx"))
    ret
}
