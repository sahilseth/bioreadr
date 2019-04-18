
get_categ <- function (mut, context = "context",
                       aAllele = "alt_allele",
                       rAllele = "ref_allele",
                       mutchar = "N"){
  
  # for null or non SNV mutation, use categ7
  if (any(c(nchar(mut[grep("alt_allele", names(mut))]),
            nchar(mut[grep("ref_allele", names(mut))])) > 1) || any(mut[grep("type", names(mut))] != "S")) {
    return(7)
  }
  
  context <- unlist(strsplit(as.character(mut["context"]), "*"))
  
  # if ref is G,
  if (mut[grep("ref_allele", names(mut))] == "G") {
    if (context[grep(mutchar, context) - 1] == "C") {
      if (isTransition(mut[grep("ref_allele", names(mut))],
                       mut[grep("alt_allele", names(mut))])) {
        return(1)
      }
      else {
        return(2)
      }
    }
  }
  
  if (mut[grep("ref_allele", names(mut))] %in% c("C", "G")) {
    if (isTransition(mut[grep("ref_allele", names(mut))],
                     mut[grep("alt_allele", names(mut))])) {
      return(3)
    }
    return(4)
  }
  if (mut[grep("ref_allele", names(mut))] %in% c("A", "T")) {
    if (isTransition(mut[grep("ref_allele", names(mut))],
                     mut[grep("alt_allele", names(mut))])) {
      return(5)
    }
    return(6)
  }
}


getEffect <- function (mut,
                       func = "func",
                       silent = c("", NA, "NA", NULL, 
                                  "unknown", "synonymous SNV", 
                                  "splicing synonymous SNV"),
                       nonsilent = c("nonsynonymous SNV", 
                                     "stopgain", "stoploss"),
                       noncoding = "."){
  
  if(mut[func] %in% silent){
    return('silent')
    
  }else if(mut[func] %in% nonsilent){
    return('nonsilent')
    
  }else if(mut[func] %in% noncoding){
    return('noncoding')
  }else{
    return('null')
  }
}

#### Called by mutect2MAF
getTumorRef <- function(mut, what = c(1, 2), trCount = "t_ref_count", 
                        rAllele = "ref_allele", aAllele = "alt_allele"){
  if(what == 1){
    if(as.numeric(mut[trCount]) > 0){
      return(mut[grep(rAllele, names(mut))])
    }else{
      return(mut[grep(aAllele, names(mut))])
    }
  }else{
    return(mut[grep(aAllele, names(mut))])
  }
}


tsv2maf <- function(x,
                    gene = "gene.refgene",
                    entrez_gene_id = "entrez_gene_id",
                    
                    chrom = "chr",
                    start_pos = "start",
                    end_pos = "end",
                    
                    func = "exonic.func",
                    ref_allele = "ref_allele",
                    alt_allele = "alt_allele",
                    
                    # AF
                    t_alt_count = "t_alt_count",
                    t_ref_count = "t_ref_count",
                    
                    sample_name = "sample_name",
                    ref_name = "sample_name",
                    sample_bam = "sample_bam",
                    
                    # mut qual
                    mutation_score = "score",
                    
                    # AA change
                    protein_change = "aachange"
                    
                    # TCGA/clinVAR etc
                    
                    
                    ){
  
  x = as.data.frame(x, stringsAsFactors = FALSE)
  
  message("seq alleles")
  Tumor_Seq_Allele1 = apply(x, 1, getTumorRef,
                            what = 1,
                            trCount = t_ref_count,
                            rAllele = ref_allele,
                            aAllele = alt_allele)
  Tumor_Seq_Allele2 = apply(x, 1,  getTumorRef, what = 2,
                            trCount = t_ref_count,
                            rAllele = ref_allele,
                            aAllele = alt_allele)
  
  maf <- data.frame(Hugo_Symbol = x[, gene],
                    #Entrez_Gene_Id = x[, entrez_gene_id],
                    Center = "IACS-MDACC",
                    NCBI_Build = "hg19",
                    
                    # position
                    Chromosome = x[, chrom],
                    Start_position = x[, start_pos],
                    End_position = x[, end_pos],
                    Strand = "+",
                    
                    # class
                    Variant_Classification = x[, func],
                    Variant_Type = "SNP",
                    
                    Reference_Allele = x[, ref_allele],
                    Tumor_Seq_Allele1 = Tumor_Seq_Allele1,
                    Tumor_Seq_Allele2 = Tumor_Seq_Allele2,
                    Match_Norm_Seq_Allele1 = "",
                    Match_Norm_Seq_Allele2 = "",
                    Tumor_Validation_Allele1 =  "",
                    Tumor_Validation_Allele2 = "",
                    
                    
                    # tumor allele counts
                    t_alt_count = x[, t_alt_count],
                    t_ref_count = x[, t_ref_count],
                    # AF:
                    t_vaf = x[, t_alt_count]/(x[, t_alt_count] + x[, t_ref_count]),
                    
                    # normal allele counts
                    
                    dbSNP_RS = x[, grep("dbsnp129", colnames(x))],
                    dbSNP_Val_Status = "bySubmitter",
                    Tumor_Sample_Barcode = x[, sample_name],
                    Matched_Norm_Sample_Barcode = x[, ref_name],
                    
                    Verification_Status = "",
                    Validation_Status = "",
                    Mutation_Status = "Somatic",
                    Sequencing_Phase = "Phase_I",
                    Sequence_Source = "Capture",
                    Validation_Method = "",
                    Score = x[, mutation_score],
                    BAM_File = x[, sample_bam],
                    Sequencer = "Illumina HiSeq",
                    effect = apply(x, 1, getEffect, func= func),
                    
                    Protein_Change = x[, protein_change],
                    
                    
                    stringsAsFactors = FALSE)
  #categ = apply(x, 1, getCateg, context=context, aAllele=alt_allele, rAllele=ref_allele))
  return(maf)
  
}


mutect_ann_to_maf = tsv2maf
