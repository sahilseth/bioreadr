
.get_value_type <- function(type){
  as_n = function(...){
    as.numeric(unlist(...))  
  }
  as_i = function(...){
    as.integer(unlist(...))  
  }
  as_c = function(...){
    as.character(unlist(...))  
  }
  
  if(type == "Float")
    return(as_n)
  if(type == "Integer")
    return(as_i)
  if(type == "String")
    return(as_c)
  
}



extract_geno <- function(vcf, tumor_name, normal_name, debug = T){
  
  num_samples = length(vcf$Samples)
  
  # geno_vars = c("GT", "AD", "AF", "DP", "GQ", "PGT", "PID", "PL", "PS")
  lst_geno = geno(vcf)#[geno_vars]
  # each of them
  # print(sapply(lst_geno, dim))
  
  # select columns which correspond to number of samples:
  # geno_vars = which(sapply(lst_geno, ncol) == num_samples) %>% names()
  geno_vars = names(lst_geno)
  rlogging::message("working on extacting:\n",
                    paste0(geno_vars, collapse = "\t"))
  # select these vars
  # lst_geno = lst_geno[geno_vars]
  
  
  # for every variable:
  # for every sample
  # unnest list into a df
  lst_geno = lapply(geno_vars, function(var){
    value_type = geno(header(vcf))[var, "Type"]
    value_number = geno(header(vcf))[var, "Number"]
    value_func = .get_value_type(value_type)
    rlogging::message(var, " ", value_type, " ", value_number, " ")
    
    df_var = geno(vcf)[[var]] %>% data.frame(check.names = F)
    # LOOP for every sample:
    # lapply(1:num_samples){
    #   
    #   
    # }
    # 
    
    # head(df_var)
    # sapply(df_var, class)
    
    
    
    if(num_samples > 1){
      # change the value
      colnames(df_var) = paste0(colnames(df_var), "_", tolower(var))
      # print(head(df_var))
      # if we have a single value, we can force typing
      # not sure what A means!!!!
      print(head(df_var))
      print(sapply(df_var, class))
      # we expect 1 number per sample, or if value_number is A
      # for the type of column
      if(value_number == 1 | value_number == "A"){
        message("     forcing type")
        df_var = apply(df_var, 2, value_func)
      }else{
        message("     skipping value-type (need to add func for this)")
      }
    }
    else{ 
      df_var = apply(df_var, 2, value_func)
      colnames(df_var) = tolower(var)
    }
    
    df_var
  })
  df_geno = do.call(cbind, lst_geno)
  # head(df_geno)
  
  # cleanup names
  if(num_samples > 1)
    colnames(df_geno) %<>% gsub(tumor_name, "t", .) %>% gsub(normal_name, "n", .)
  df_geno
  
}


extract_fr_single <- function(vcf, samplename){
  
  sample_name = colnames(geno(vcf)[[1]])
  unnest_fr <- function(x, suffix){
    tmp1 = x[,1] %>% do.call(rbind, .)
    # add ref and alt
    colnames(tmp1) = c("ref", "alt")
    # then add the suffix
    colnames(tmp1) = paste0(colnames(tmp1), "_", suffix)
    tmp1
  }
  
  df_f1r2 = unnest_fr(geno(vcf)$F1R2, "f1r2")
  df_f2r1 = unnest_fr(geno(vcf)$F2R1, "f2r1")
  df_fr = cbind(df_f1r2, df_f2r1)
  print(colnames(df_fr))
  df_fr %>% wranglr::to_df() %>%
    mutate(af_f1r2 = alt_f1r2/(ref_f1r2 + alt_f1r2),
           af_f2r1 = alt_f2r1/(ref_f2r1 + alt_f2r1))
}

extract_fr_pair <- function(vcf, tumor_name, normal_name){
  
  check_args()
  unnest_fr <- function(x, suffix){
    tmp1 = x[,1] %>% do.call(rbind, .)
    tmp2 = x[,2] %>% do.call(rbind, .)
    tmp3 = cbind(tmp1, tmp2)
    colnames(tmp3) = c(paste0(colnames(x)[1], c("_ref", "_alt")),
                       paste0(colnames(x)[2], c("_ref", "_alt")))
    colnames(tmp3) = paste0(colnames(tmp3), "_", suffix)
    tmp3
  }
  
  df_f1r2 = unnest_fr(geno(vcf)$F1R2, "f1r2")
  df_f2r1 = unnest_fr(geno(vcf)$F2R1, "f2r1")
  df_fr = cbind(df_f1r2, df_f2r1)
  colnames(df_fr) %<>% gsub(tumor_name, "t", .) %>% gsub(normal_name, "n", .)
  df_fr %>% wranglr::to_df() %>%
    mutate(t_af_f1r2 = t_alt_f1r2/(t_ref_f1r2 + t_alt_f1r2),
           t_af_f2r1 = t_alt_f2r1/(t_ref_f2r1 + t_alt_f2r1))
}


uniquefy_coloumn_names <- function(lst_dfs, names_dfs){
  
  df_col_names = lapply(1:length(lst_dfs), function(i){
    df = lst_dfs[[i]]
    cols = colnames(df)
    data.frame(column_name = cols, df_name = names_dfs[i], stringsAsFactors = F)
  }) %>% bind_rows()
  df_col_names %<>% group_by(column_name) %>%
    mutate(num = n()) %>% ungroup() %>%
    mutate(column_name2 = ifelse(num > 1, paste0(column_name, "_", df_name), column_name))
}

# http://adv-r.had.co.nz/S3.html
to_df <- function (x, ...) {
  UseMethod("to_df", x)
}

to_df.CollapsedVCF <- function(vcf, tumor_name, normal_name){
  
  rlogging::message("input is of class ", class(vcf))
  
  rlogging::message("extracing rowranges")
  # df_fixed = vcf@fixed %>% as.data.frame() %>% clean_names()
  df_fixed = rowRanges(vcf) %>% data.frame() %>% clean_names()
  
  rlogging::message("extracing info")
  df_info = info(vcf) %>% as.data.frame() %>% clean_names()
  
  num_samples = length(vcf$Samples)
  rlogging::message("extracing geno")
  # dont include SB (that has a lot more columns)
  # geno_vars = c("GT", "AD", "AF", "DP", "GQ", "PGT", "PID", "PL", "PS")
  if(num_samples > 1)
    df_geno = extract_geno(vcf, tumor_name, normal_name)
  else
    df_geno = extract_geno(vcf)
  
  
  # rlogging::message("extracing fr")
  # if(num_samples > 1)
  #   df_fr = extract_fr_pair(vcf, tumor_name, normal_name)
  # else
  #   df_fr = extract_fr_single(vcf)
  # # still need a way to sort TN stuff
  # # df_vcf = cbind(df_fixed, df_info)
  
  rlogging::message("merging dfs (unique column names)")
  df_coloumn_names = uniquefy_coloumn_names(
    list(df_fixed, df_info, df_geno, df_fr), 
    c("fixed", "info", "geno", "fr"))
  
  rlogging::message("cbind")
  df_vcf = cbind(df_fixed, df_info, df_geno, df_fr)
  colnames(df_vcf) = df_coloumn_names$column_name2
  # df_vcf
  
  rlogging::message("cleanup column names")
  df_vcf$start = as.integer(df_vcf$start)
  # df_vcf$start = as.integer(df_vcf$start)
  # df_vcf$af = as.numeric(df_vcf$af)
  # df_vcf$seqnames = as.character(df_vcf$seqnames) %>% paste0("chr", .)
  # df_vcf$alt = df_vcf$alt %>% sapply(as.character)
  
  df_vcf
}


# example -------
if(FALSE){
  library(pacman)
  p_load(VariantAnnotation, dplyr, magrittr, readr, janitor)
  p_load(flowr)
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/mutect2_vcf_qc.R')
  fl_m2 = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/ssm/m1_m2_ir/WEX-334187-T___matched_combvcf.vcf.gz"
  vcf_m2_f2 = VariantAnnotation::readVcf(fl_m2, "hg19")
  
  vcf_m2_f2_snv = vcf_m2_f2[isSNV(vcf_m2_f2), ]
  vcf = vcf_m2_f2_snv
  geno(header(vcf))
  # debug(to_df.CollapsedVCF)
  df_m2_f2 = to_df(vcf_m2_f2_snv, 
                   tumor_name = '334187-T',
                   normal_name = '334187-N')
  
  
}

