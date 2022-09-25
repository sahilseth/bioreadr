


#' Title
#'
#' @param x 
#'
#' @return
#' @export
#' 
#' @details 
#' The pVACseq pipeline will write its results in separate folders depending on which prediction algorithms were chosen:
#' MHC_Class_I: for MHC class I prediction algorithms
#' MHC_Class_II: for MHC class II prediction algorithms
#' combined: If both MHC class I and MHC class II prediction algorithms were run, this folder combines the neoeptiope predictions from both
#' 
#' Each folder will contain the same list of output files (listed in the order created):
#'   File Name Description <sample_name>.tsv
#'   An intermediate file with variant, transcript, coverage, vaf, and expression information parsed from the input files.
#'   <sample_name>.tsv_<chunks> (multiple
#'   The above file but split into smaller chunks for easier processing with IEDB.
#'   
#'   <sample_name>.all_epitopes.tsv
#'   A list of all predicted epitopes and their binding affinity scores, with additional variant information from the <sample_name>.tsv.
#'   <sample_name>.filtered.tsv 
#'   The above file after applying all filters, with cleavage site and stability predictions added.
#'   <sample_name>.filtered.condensed.ranked.tsv
#'   A condensed version of the filtered TSV with only the most important columns remaining, with a priority score for each neoepitope candidate added.

# this one uses fread instead
read_pvacseq2 <- function(oprefix,
                         seq_sample_id, 
                         path_patient_id){
  p_load(janitor)

  message(" all ", appendLF = F)
  glue("{oprefix}.all_epitopes.tsv") %>% file.exists()
  # protein_position: to take care of indels
  df_all = glue("{oprefix}.all_epitopes.tsv") %>% 
    data.table::fread(data.table = F) %>% 
    dplyr::mutate(seq_sample_id = seq_sample_id, 
                  path_patient_id = path_patient_id) %>% 
    clean_names() %>% 
    dplyr::mutate(protein_position = as.character(protein_position),
           chromosome = as.character(chromosome),
           reference = as.character(reference))
  
  message("filtered ", appendLF = F)
  df_filt = glue("{oprefix}.filtered.tsv") %>%
    data.table::fread(data.table = F) %>%
      dplyr::mutate(seq_sample_id = seq_sample_id,
                    path_patient_id = path_patient_id) %>% 
      clean_names() %>%
        # to take care of indels
      dplyr::mutate(
        protein_position = as.character(protein_position),
        chromosome = as.character(chromosome),
        reference = as.character(reference),
        transcript = as.character(transcript),
        ensembl_gene_id = as.character(ensembl_gene_id),
        variant = as.character(variant),
        mutation = as.character(mutation),
        gene_name = as.character(gene_name),
        hla_allele = as.character(hla_allele),
        mt_epitope_seq = as.character(mt_epitope_seq),
        wt_epitope_seq = as.character(wt_epitope_seq),
        best_mt_score_method = as.character(best_mt_score_method),
        variant_type = as.character(variant_type))

  message("filtered ranked")
  df_filt_rnk = glue("{oprefix}.filtered.condensed.ranked.tsv") %>%
    data.table::fread(data.table = F) %>%
      dplyr::mutate(
        seq_sample_id = seq_sample_id,
        path_patient_id = path_patient_id
      ) %>% clean_names() %>% 
    # to take care of indels
    dplyr::mutate(
      protein_position = as.character(protein_position),
      gene_name = as.character(gene_name),
      hla_allele = as.character(hla_allele),
      mt_epitope_seq = as.character(mt_epitope_seq),
      mutation = as.character(mutation))
  
  
  list(df_all = df_all, 
       df_filt = df_filt,
       df_filt_rnk = df_filt_rnk)
}

read_pvacseq <- function(oprefix,
                         seq_sample_id, 
                         path_patient_id){
  
  df_all = read_tsv(glue("{oprefix}.all_epitopes.tsv"), 
                    col_types = cols(
                      Chromosome = col_character(),
                      Start = col_double(),
                      Stop = col_double(),
                      Reference = col_character(),
                      Variant = col_character(),
                      Transcript = col_character(),
                      `Transcript Support Level` = col_logical(),
                      `Ensembl Gene ID` = col_character(),
                      `Variant Type` = col_character(),
                      Mutation = col_character(),
                      `Protein Position` = col_double(),
                      `Gene Name` = col_character(),
                      HGVSc = col_logical(),
                      HGVSp = col_logical(),
                      `HLA Allele` = col_character(),
                      `Peptide Length` = col_double(),
                      `Sub-peptide Position` = col_double(),
                      `Mutation Position` = col_double(),
                      `MT Epitope Seq` = col_character(),
                      `WT Epitope Seq` = col_character(),
                      `Best MT Score Method` = col_character(),
                      `Best MT Score` = col_double(),
                      `Corresponding WT Score` = col_double(),
                      `Corresponding Fold Change` = col_double(),
                      `Tumor DNA Depth` = col_double(),
                      `Tumor DNA VAF` = col_double(),
                      `Tumor RNA Depth` = col_logical(),
                      `Tumor RNA VAF` = col_logical(),
                      `Normal Depth` = col_logical(),
                      `Normal VAF` = col_logical(),
                      `Gene Expression` = col_logical(),
                      `Transcript Expression` = col_logical(),
                      `Median MT Score` = col_double(),
                      `Median WT Score` = col_double(),
                      `Median Fold Change` = col_double(),
                      `MHCnuggetsI WT Score` = col_double(),
                      `MHCnuggetsI MT Score` = col_double(),
                      cterm_7mer_gravy_score = col_double(),
                      max_7mer_gravy_score = col_double(),
                      difficult_n_terminal_residue = col_logical(),
                      c_terminal_cysteine = col_logical(),
                      c_terminal_proline = col_logical(),
                      cysteine_count = col_double(),
                      n_terminal_asparagine = col_logical(),
                      asparagine_proline_bond_count = col_double()
                    )) %>% 
    dplyr::mutate(seq_sample_id = seq_sample_id, 
                  path_patient_id = path_patient_id)
  
  df_filt = read_tsv(glue("{oprefix}.filtered.tsv"), 
                     col_types = cols(
                       Chromosome = col_character(),
                       Start = col_double(),
                       Stop = col_double(),
                       Reference = col_character(),
                       Variant = col_character(),
                       Transcript = col_character(),
                       `Transcript Support Level` = col_logical(),
                       `Ensembl Gene ID` = col_character(),
                       `Variant Type` = col_character(),
                       Mutation = col_character(),
                       `Protein Position` = col_double(),
                       `Gene Name` = col_character(),
                       HGVSc = col_logical(),
                       HGVSp = col_logical(),
                       `HLA Allele` = col_character(),
                       `Peptide Length` = col_double(),
                       `Sub-peptide Position` = col_double(),
                       `Mutation Position` = col_double(),
                       `MT Epitope Seq` = col_character(),
                       `WT Epitope Seq` = col_character(),
                       `Best MT Score Method` = col_character(),
                       `Best MT Score` = col_double(),
                       `Corresponding WT Score` = col_double(),
                       `Corresponding Fold Change` = col_double(),
                       `Tumor DNA Depth` = col_double(),
                       `Tumor DNA VAF` = col_double(),
                       `Tumor RNA Depth` = col_logical(),
                       `Tumor RNA VAF` = col_logical(),
                       `Normal Depth` = col_logical(),
                       `Normal VAF` = col_logical(),
                       `Gene Expression` = col_logical(),
                       `Transcript Expression` = col_logical(),
                       `Median MT Score` = col_double(),
                       `Median WT Score` = col_double(),
                       `Median Fold Change` = col_double(),
                       `MHCnuggetsI WT Score` = col_double(),
                       `MHCnuggetsI MT Score` = col_double(),
                       cterm_7mer_gravy_score = col_double(),
                       max_7mer_gravy_score = col_double(),
                       difficult_n_terminal_residue = col_logical(),
                       c_terminal_cysteine = col_logical(),
                       c_terminal_proline = col_logical(),
                       cysteine_count = col_double(),
                       n_terminal_asparagine = col_logical(),
                       asparagine_proline_bond_count = col_double()
                     )) %>% 
    dplyr::mutate(seq_sample_id = seq_sample_id, 
                  path_patient_id = path_patient_id)
  
  df_filt_rnk = read_tsv(glue("{oprefix}.filtered.condensed.ranked.tsv"), 
                         col_types = cols(
                           `Gene Name` = col_character(),
                           Mutation = col_character(),
                           `Protein Position` = col_double(),
                           HGVSc = col_logical(),
                           HGVSp = col_logical(),
                           `HLA Allele` = col_character(),
                           `Mutation Position` = col_double(),
                           `MT Epitope Seq` = col_character(),
                           `Median MT Score` = col_double(),
                           `Median WT Score` = col_double(),
                           `Median Fold Change` = col_double(),
                           `Best MT Score` = col_double(),
                           `Corresponding WT Score` = col_double(),
                           `Corresponding Fold Change` = col_double(),
                           `Tumor DNA Depth` = col_double(),
                           `Tumor DNA VAF` = col_double(),
                           `Tumor RNA Depth` = col_logical(),
                           `Tumor RNA VAF` = col_logical(),
                           `Gene Expression` = col_logical(),
                           Rank = col_double()
                         )) %>% 
    dplyr::mutate(seq_sample_id = seq_sample_id, 
                  path_patient_id = path_patient_id)
  
  list(df_all = df_all, 
       df_filt = df_filt,
       df_filt_rnk = df_filt_rnk)
}



# trim HLA to 4 digits
# x = "DQB1*04:02:01"
# x = "DQA1*02:01"
# x = "B*44:02:01"
hla_trim_d4 <- function(x){
  lapply(x, function(xi){
    xs = strsplit(xi, split = ":")[[1]]
    # get first two elements
    paste0(xs[1:2], collapse = ":")
  })
    
}












# END
