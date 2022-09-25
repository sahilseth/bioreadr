

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
