# recode states ---------
# total/minor CN encoding -------
# https://github.com/mskcc/facets/issues/118

# https://github.com/mskcc/facets/issues/62
# Filtering and interpreting results
# tcn.em	lcn.em	type
# 0	0	DEL
# 0	NA	DEL
# 1	0	HETD
# 1	NA	HETD
# 2	0	LOH
# 2	1	NEUTR
# 2	NA	NEUTR/Unknown
# 2+	0	DUP-LOH
# 2+	1+	DUP
# 2+	NA	DUP-[LOH?]
recode_cn.df_seg <- function(df_seg){
  wranglr::expect_columns(df_seg, c("total_cn", "minor_cn"))
  df_seg = mutate(df_seg,
                  cna_status = case_when(
                    total_cn == 0 & minor_cn == 0 ~ "DEL",
                    total_cn == 0 & is.na(minor_cn) ~ "DEL",
                    
                    total_cn == 1 & minor_cn == 0 ~ "HETD",
                    total_cn == 1 & is.na(minor_cn) == 0 ~ "HETD",
                    
                    total_cn == 2 & minor_cn == 0 ~ "LOH",
                    
                    total_cn == 2 & minor_cn == 1 ~ "NEU",
                    total_cn == 2 & is.na(minor_cn) ~ "NEU/UNK",
                    
                    total_cn > 2 & total_cn < 4 & minor_cn == 0 ~ "DUP-LOH",
                    total_cn > 2 & total_cn < 4 & is.na(minor_cn) ~ "DUP-LOH",
                    
                    total_cn > 2 & minor_cn == 1 ~ "DUP",
                    total_cn > 2 & minor_cn >= 1 ~ "GAIN",
                    
                    # there are several cases where # of het sites are not available
                    # therefore
                    total_cn >= 4 ~ "AMP"
                  ),
                  cna_status = factor(cna_status, levels = c("DEL", "HETD",
                                                             "LOH", "NEU",
                                                             "NEU/UNK", "DUP-LOH",
                                                             "DUP",
                                                             "GAIN", "AMP")))
  
  # attr(df_seg$cna_status, "tcn") = df_seg$total_cn
  df_seg
}
recode_cn.facets = recode_cn.df_seg
# favoring deletions to resolve ties
cna_state_order = c("DEL" = -2.1, "HETD" = -1.1, "LOH" = -0.5, 
                    "NEU" = 0, "NEU/UNK" = 0, 
                    "DUP-LOH" = 0.5, "DUP" = 1, "GAIN" = 2, "AMP" = 3)


# resolve status in case of multiple per gene
# gets the max
resolve_cna_status <- function(x){
  if (length(x) == 1) {
    return(x)
  }
  # x = df_igv_seg$cna_status[1:3]
  xs = abs(cna_state_order[x])
  # if unresolves, just take one
  x[which(xs == max(xs))][1]
}
# resolve_cna_status(c("DEL", "GAIN"))
# resolve_cna_status(c("DEL", "AMP"))



colors_cna_status.facets <- function(){
  c("DEL" = "#762a83", #dark purple
    "HEMIZYGOTE" = "#af8dc3", #light purple
    "LOH" = "seagreen", #light purple
    "NEU/UNK" = "gray88",
    "NEU" = "gray88",
    "DUP-LOH" = "#7fbf7b", #light green
    "GAIN" = "#1b7837", #dark green
    "AMP" = "red4",
    '(Missing)' = "white",
    "'(Missing)'" = "white")
}
