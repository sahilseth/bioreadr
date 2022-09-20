# https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format
# The process for modifying a protected MAF into a somatic MAF is as follows:
#   
#   Aliquot Selection: only one tumor-normal pair are selected for each tumor sample based on the plate number, sample type, analyte type and other features extracted from tumor TCGA aliquot barcode.
# Low quality variant filtering and germline masking:
#   Variants with Mutation_Status != 'Somatic' or GDC_FILTER = 'Gapfiller', 'ContEst', 'multiallelic', 'nonselectedaliquot', 'BCR_Duplicate' or 'BadSeq' are removed.
# Remaining variants with GDC_Valid_Somatic = True are included in the Somatic MAF.
# Remaining variants with FILTER != 'panel_of_normals' or PASS are removed. Note that the FILTER != panel_of_normals value is only relevant for the variants generated from the MuTect2 pipeline.
# Remaining variants with MC3_Overlap = True are included in the Somatic MAF.
# Remaining variants with GDC_FILTER = 'ndp', 'NonExonic', 'bitgt', 'gdc_pon' are removed.
# Remaining variants with SOMATIC != null are included in the Somatic MAF.
# Remaining variants with dbSNP_RS = 'novel' or null are included in the Somatic MAF.
# Remaining variants are removed.
# Removal of the following columns:
#   vcf_region
# vcf_info
# vcf_format
# vcf_tumor_gt
# vcf_normal_gt
# GDC_Valid_Somatic
# Set values to be blank in the following columns that may contain information about germline genotypes:
# Match_Norm_Seq_Allele1
# Match_Norm_Seq_Allele2
# Match_Norm_Validation_Allele1
# Match_Norm_Validation_Allele2
# n_ref_count
# n_alt_count


# Protected MAF File Structure 
# The table below describes the columns in a protected MAF and their definitions. Note that the somatic (open-access) MAF structure is the same except for having the last six columns removed.
# 
# Column	Description
# 1 - Hugo_Symbol	  HUGO symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene
# 2 - Entrez_Gene_Id	  Entrez gene ID (an integer). "0" is used for regions that do not correspond to a gene region or Ensembl ID
# 3 - Center	One or more genome sequencing center reporting the variant
# 4 - NCBI_Build	The reference genome used for the alignment (GRCh38)
# 5 - Chromosome	The affected chromosome (chr1)
# 6 - Start_Position	Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate
# 7 - End_Position	Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate
# 8 - Strand	Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+'
# 9 - Variant_Classification	Translational effect of variant allele
# 10 - Variant_Type	Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)
# 11 - Reference_Allele	The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion
# 12 - Tumor_Seq_Allele1	Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
# 13 - Tumor_Seq_Allele2	Tumor sequencing (discovery) allele 2
# 14 - dbSNP_RS	The rs-IDs from the   dbSNP database, "novel" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases
# 15 - dbSNP_Val_Status	The dbSNP validation status is reported as a semicolon-separated list of statuses. The union of all rs-IDs is taken when there are multiple
# 16 - Tumor_Sample_Barcode	Aliquot barcode for the tumor sample
# 17 - Matched_Norm_Sample_Barcode	Aliquot barcode for the matched normal sample
# 18 - Match_Norm_Seq_Allele1	Primary data genotype. Matched normal sequencing allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)
# 19 - Match_Norm_Seq_Allele2	Matched normal sequencing allele 2
# 20 - Tumor_Validation_Allele1	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
# 21 - Tumor_Validation_Allele2	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2
# 22 - Match_Norm_Validation_Allele1	Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)
# 23 - Match_Norm_Validation_Allele2	Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2 (cleared in somatic MAF)
# 24 - Verification_Status	Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing
# 25 - Validation_Status	Second pass results from orthogonal technology
# 26 - Mutation_Status	An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field
# 27 - Sequencing_Phase	TCGA sequencing phase (if applicable). Phase should change under any circumstance that the targets under consideration change
# 28 - Sequence_Source	Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the SRA 1.5 library_strategy field values. This subset matches those used at CGHub
# 29 - Validation_Method	The assay platforms used for the validation call
# 30 - Score	Not in use
# 31 - BAM_File	Not in use
# 32 - Sequencer	Instrument used to produce primary sequence data
# 33 - Tumor_Sample_UUID	GDC aliquot UUID for tumor sample
# 34 - Matched_Norm_Sample_UUID	GDC aliquot UUID for matched normal sample
# 35 - HGVSc	The coding sequence of the variant in HGVS recommended format
# 36 - HGVSp	The protein sequence of the variant in HGVS recommended format. "p.=" signifies no change in the protein
# 37 - HGVSp_Short	Same as the HGVSp column, but using 1-letter amino-acid codes
# 38 - Transcript_ID	  Ensembl ID of the transcript affected by the variant
# 39 - Exon_Number	The exon number (out of total number)
# 40 - t_depth	Read depth across this locus in tumor BAM
# 41 - t_ref_count	Read depth supporting the reference allele in tumor BAM
# 42 - t_alt_count	Read depth supporting the variant allele in tumor BAM
# 43 - n_depth	Read depth across this locus in normal BAM
# 44 - n_ref_count	Read depth supporting the reference allele in normal BAM (cleared in somatic MAF)
# 45 - n_alt_count	Read depth supporting the variant allele in normal BAM (cleared in somatic MAF)
# 46 - all_effects	A semicolon delimited list of all possible variant effects, sorted by priority ([Symbol,Consequence,HGVSp_Short,Transcript_ID,RefSeq,HGVSc,Impact,Canonical,Sift,PolyPhen,Strand])
# 47 - Allele	The variant allele used to calculate the consequence
# 48 - Gene	Stable Ensembl ID of affected gene
# 49 - Feature	Stable Ensembl ID of feature (transcript, regulatory, motif)
# 50 - Feature_type	Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (or blank)
# 51 - One_Consequence	The single consequence of the canonical transcript in   sequence ontology terms
# 52 - Consequence	Consequence type of this variant;   sequence ontology terms
# 53 - cDNA_position	Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA
# 54 - CDS_position	Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
# 55 - Protein_position	Relative position of affected amino acid in protein. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
# 56 - Amino_acids	Only given if the variation affects the protein-coding sequence
# 57 - Codons	The alternative codons with the variant base in upper case
# 58 - Existing_variation	Known identifier of existing variation
# 59 - ALLELE_NUM	Allele number from input; 0 is reference, 1 is first alternate etc.
# 60 - DISTANCE	Shortest distance from the variant to transcript
# 61 - TRANSCRIPT_STRAND	The DNA strand (1 or -1) on which the transcript/feature lies
# 62 - SYMBOL	The gene symbol
# 63 - SYMBOL_SOURCE	The source of the gene symbol
# 64 - HGNC_ID	Gene identifier from the HUGO Gene Nomenclature Committee if applicable
# 65 - BIOTYPE	Biotype of transcript
# 66 - CANONICAL	A flag (YES) indicating that the VEP-based canonical transcript, the longest translation, was used for this gene. If not, the value is null
# 67 - CCDS	The   CCDS identifier for this transcript, where applicable
# 68 - ENSP	The Ensembl protein identifier of the affected transcript
# 69 - SWISSPROT	  UniProtKB/Swiss-Prot accession
# 70 - TREMBL	UniProtKB/TrEMBL identifier of protein product
# 71 - UNIPARC	UniParc identifier of protein product
# 72 - RefSeq	RefSeq identifier for this transcript
# 73 - SIFT	The   SIFT prediction and/or score, with both given as prediction (score)
# 74 - PolyPhen	The   PolyPhen prediction and/or score
# 75 - EXON	The exon number (out of total number)
# 76 - INTRON	The intron number (out of total number)
# 77 - DOMAINS	The source and identifier of any overlapping protein domains
# 78 - GMAF	Non-reference allele and frequency of existing variant in   1000 Genomes
# 79 - AFR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined African population
# 80 - AMR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined American population
# 81 - ASN_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population
# 82 - EAS_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population
# 83 - EUR_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined European population
# 84 - SAS_MAF	Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population
# 85 - AA_MAF	Non-reference allele and frequency of existing variant in   NHLBI-ESP African American population
# 86 - EA_MAF	Non-reference allele and frequency of existing variant in NHLBI-ESP European American population
# 87 - CLIN_SIG	Clinical significance of variant from dbSNP
# 88 - SOMATIC	Somatic status of each ID reported under Existing_variation (0, 1, or null)
# 89 - PUBMED	Pubmed ID(s) of publications that cite existing variant
# 90 - MOTIF_NAME	The source and identifier of a transcription factor binding profile aligned at this position
# 91 - MOTIF_POS	The relative position of the variation in the aligned TFBP
# 92 - HIGH_INF_POS	A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (Y, N, or null)
# 93 - MOTIF_SCORE_CHANGE	The difference in motif score of the reference and variant sequences for the TFBP
# 94 - IMPACT	The impact modifier for the consequence type
# 95 - PICK	Indicates if this block of consequence data was picked by VEP's   pick feature (1 or null)
# 96 - VARIANT_CLASS	Sequence Ontology variant class
# 97 - TSL	  Transcript support level, which is based on independent RNA analyses
# 98 - HGVS_OFFSET	Indicates by how many bases the HGVS notations for this variant have been shifted
# 99 - PHENO	Indicates if existing variant is associated with a phenotype, disease or trait (0, 1, or null)
# 100 - MINIMISED	Alleles in this variant have been converted to minimal representation before consequence calculation (1 or null)
# 101 - ExAC_AF	Global Allele Frequency from   ExAC
# 102 - ExAC_AF_Adj	Adjusted Global Allele Frequency from ExAC
# 103 - ExAC_AF_AFR	African/African American Allele Frequency from ExAC
# 104 - ExAC_AF_AMR	American Allele Frequency from ExAC
# 105 - ExAC_AF_EAS	East Asian Allele Frequency from ExAC
# 106 - ExAC_AF_FIN	Finnish Allele Frequency from ExAC
# 107 - ExAC_AF_NFE	Non-Finnish European Allele Frequency from ExAC
# 108 - ExAC_AF_OTH	Other Allele Frequency from ExAC
# 109 - ExAC_AF_SAS	South Asian Allele Frequency from ExAC
# 110 - GENE_PHENO	Indicates if gene that the variant maps to is associated with a phenotype, disease or trait (0, 1, or null)
# 111 - FILTER	Copied from input VCF. This includes filters implemented directly by the variant caller and other external software used in the DNA-Seq pipeline. See below for additional details.
# 112 - CONTEXT	The reference allele per VCF specs, and its five flanking base pairs
# 113 - src_vcf_id	GDC UUID for the input VCF file
# 114 - tumor_bam_uuid	GDC UUID for the tumor bam file
# 115 - normal_bam_uuid	GDC UUID for the normal bam file
# 116 - case_id	GDC UUID for the case
# 117 - GDC_FILTER	GDC filters applied universally across all MAFs
# 118 - COSMIC	Overlapping COSMIC variants
# 119 - MC3_Overlap	Indicates whether this region overlaps with an MC3 variant for the same sample pair
# 120 - GDC_Validation_Status	GDC implementation of validation checks. See notes section (#5) below for details
# 121 - GDC_Valid_Somatic	True or False (not in somatic MAF)
# 122 - vcf_region	Colon separated string containing the CHROM, POS, ID, REF, and ALT columns from the VCF file (e.g., chrZ:20:rs1234:A:T) (not in somatic MAF)
# 123 - vcf_info	INFO column from VCF (not in somatic MAF)
# 124 - vcf_format	FORMAT column from VCF (not in somatic MAF)
# 125 - vcf_tumor_gt	Tumor sample genotype column from VCF (not in somatic MAF)
# 126 - vcf_normal_gt	Normal sample genotype column from VCF (not in somatic MAF)



# 11 - Reference_Allele	The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion
# 12 - Tumor_Seq_Allele1	Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases
# 13 - Tumor_Seq_Allele2	Tumor sequencing (discovery) allele 2
#### Called by mutect2MAF





# parse VCF into a DF (parse_somatic_vcf)
# filter it
# convert to maf

# https://bioconductor.org/packages/3.7/bioc/vignettes/maftools/inst/doc/maftools.html#8_variant_annotations

##   chromsome    start      end     ref alt Tumor_Sample_Barcode
## 1      chr4 55589774 55589774       A   G               fake_1
## 2      chr4 55599321 55599321       A   T               fake_2
## 3      chr4 55599332 55599332       G   T               fake_3
## 4      chr4 55599320 55599320       G   T               fake_4
## 5     chr15 41961117 41961123 TGGCTAA   -               fake_4
## 6      chr4 55599320 55599320       G   T               fake_5



to_maf <- function (x, ...) {
  
  UseMethod("to_maf", x)
  
}

# the most basic of df
to_maf.mutect <- function(){
  
}

# the most basic of df
to_maf.mutect_annovar <- function(){
  
}

# the most basic of df
to_maf.mutect2 <- function(){
  
}

# the most basic of df
# this is dummy correct, not the default output of pcgr
to_maf.pcgr <- function(x, maftools = FALSE) {
  # x = df_mut
  x %<>% mutate(
    consequence2 = recode_variant_classification.pcgr(consequence),
    t_alt_count = round(tdp*taf, 2),
    t_ref_count = tdp - t_alt_count)
  print(table(x$consequence2))
  maf = to_maf.data.frame(x, 
                          gene = "symbol", 
                          entrez_gene_id = "entrez_id",

                          chrom = "chrom",
                          start_pos = "pos",
                          end_pos = "pos",

                          func = "consequence2",
                          ref_allele = "ref",
                          alt_allele = "alt",

                          # AF
                          t_alt_count = "t_alt_count",
                          t_ref_count = "t_ref_count",

                          sample_name = "name",
                          ref_name = "name",
                          sample_bam = "name",

                          # mut qual
                          mutation_score = "m2_qual",

                          # AA change
                          aa_change = "protein_change")

    if(maftools)
      maf = maftools::read.maf(maf)

  return(maf)
}

# the most basic of df
to_maf.data.frame <- function(x,
                              gene = "gene.knowngene",
                              entrez_gene_id = "entrez_gene_id",

                              chrom = "chr",
                              start_pos = "start",
                              end_pos = "end",

                              func = "exonicfunc.knowngene",
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
                              aa_change = "aachange"

                              # TCGA/clinVAR etc


) {
  p_load(futile.logger)

  if (is.data.frame(x)) {
    x = as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x = read_tsv(x, col_types = cols(.default = col_character())) %>%
      data.frame(stringsAsFactors = FALSE)
  }

  cols_data = colnames(x)
  cols_expected = c(
    gene, entrez_gene_id,
    chrom, start_pos, end_pos,
    func, ref_allele, alt_allele,
    t_alt_count, t_ref_count,
    sample_bam, ref_name, sample_name,
    mutation_score, aa_change
  )
  cols_missing = cols_expected[!cols_expected %in% cols_data]
  if (length(cols_missing) > 0) {
    str_cols_missing = paste0(cols_missing, collapse = "\n")
    str_cols_data <- paste0(cols_data, collapse = "\n")
    flog.warn(glue("These columns are missing: \n",
      "{str_cols_missing}\n\n",
      "Available options are:\n",
      "{str_cols_data}"))
    stop("fix missing columns")
  }
  x[, t_alt_count] = as.integer(x[, t_alt_count])
  x[, t_ref_count] = as.integer(x[, t_ref_count])


  flog.info("> variant type")
  x$.variant_type = get_variant_type(x[, ref_allele], x[, alt_allele])

  flog.info("> seq alleles tumor allele1")
  Tumor_Seq_Allele1 = apply(x, 1, getTumorRef,
    what = 1,
    trCount = t_ref_count,
    rAllele = ref_allele,
    aAllele = alt_allele,
    variant_type = ".variant_type"
  )
  flog.info("> seq alleles tumor allele2")
  Tumor_Seq_Allele2 = apply(x, 1, getTumorRef,
    what = 2,
    trCount = t_ref_count,
    rAllele = ref_allele,
    aAllele = alt_allele,
    variant_type = ".variant_type"
  )
  flog.info("> seq alleles ref allele")
  Reference_Allele = get_reference_allele(x[, ref_allele], x$.variant_type)


  flog.info("creating df")
  maf <- data.frame(
    Hugo_Symbol = x[, gene],
    Entrez_Gene_Id = x[, entrez_gene_id],
    Center = "IACS-MDACC",
    NCBI_Build = "hg19",

    # position
    Chromosome = x[, chrom],
    Start_position = x[, start_pos],
    End_position = x[, end_pos],
    Strand = "+",

    # class
    Variant_Classification = x[, func],
    Variant_Type = x$.variant_type,

    Reference_Allele = Reference_Allele,
    Tumor_Seq_Allele1 = Tumor_Seq_Allele1,
    Tumor_Seq_Allele2 = Tumor_Seq_Allele2,
    Match_Norm_Seq_Allele1 = "",
    Match_Norm_Seq_Allele2 = "",
    Tumor_Validation_Allele1 = "",
    Tumor_Validation_Allele2 = "",


    # tumor allele counts
    t_alt_count = x[, t_alt_count],
    t_ref_count = x[, t_ref_count],
    # AF:
    t_vaf = x[, t_alt_count] / (x[, t_alt_count] + x[, t_ref_count]),

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
    effect = apply(x, 1, getEffect, func = func),

    Protein_Change = x[, aa_change],


    stringsAsFactors = FALSE
  )


  # categ = apply(x, 1, getCateg, context=context, aAllele=alt_allele, rAllele=ref_allele))
  return(maf)
}
mutect_ann_to_maf = tsv2maf = to_maf.data.frame



#' to maf
#'
#' @param df is a output of parse_somatic_vcf
#'
.to_maf.vcf <- function(df){
  
  parse_somatic_vcf
  
  
}



.get_categ <- function (mut, 
                        context = "context",
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




getTumorRef <- function(mut, what = c(1, 2), 
                        trCount = "t_ref_count", 
                        rAllele = "ref_allele", 
                        aAllele = "alt_allele", 
                        variant_type){
  variant_type = mut[variant_type]
  if(what == 1){
    if(variant_type == "DEL")
      return("-")
    if(variant_type %in%  c("DNP", "TNP", "ONP"))
      return(aAllele)

    # handle SNPs
    # if tumor ref counts > 0, get ref_allele
    if(as.numeric(mut[trCount]) > 0){
      # return(mut[grep(rAllele, names(mut))])
      return(mut[rAllele])
    # if tumor ref counts == 0, there is no tumor reference, only ALT allele
    }else{
      return(mut[aAllele])
      # return(mut[grep(aAllele, names(mut))])
    }

  }else{
    # case of allele2, always return alt allele
    # fill missing allele, in case of complex rearrangements
    if(variant_type %in% c("DEL", "INS", "DNP", "TNP", "ONP"))
      return("")
    else
      return(mut[aAllele])
  }
}

# Variant_Type	Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)
get_variant_type <- function(ref, alt){
  variant_type <- case_when(
    nchar(ref) > nchar(alt) ~ "DEL",
    nchar(ref) < nchar(alt) ~ "INS",
    nchar(ref) == nchar(alt) & nchar(ref) == 1 ~ "SNP",
    nchar(ref) == nchar(alt) & nchar(ref) == 2 ~ "DNP",
    nchar(ref) == nchar(alt) & nchar(ref) == 3 ~ "TNP",
    nchar(ref) == nchar(alt) & nchar(ref) > 3 ~ "ONP", 
    TRUE ~ "ONP")
}

get_reference_allele <- function(ref, variant_type){
  # maf2vcf throws an error that ref bases do not match
  ref_del = substring(ref, 2)
  ref <- case_when(
    variant_type == "INS" ~ "-",
    variant_type == "DEL" ~ ref_del,
    TRUE ~ ref)
}






# https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/


# variant_type:
# Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. 
# ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)

# https://github.com/mskcc/vcf2maf/blob/7f3bd61799c45adb6afbe51e0f21210753904fd5/vcf2maf.pl
# my ( $effect, $var_type, $inframe ) = @_;
# return "Splice_Site" if( $effect =~ /^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)$/ );
# return "Nonsense_Mutation" if( $effect eq 'stop_gained' );
# return "Frame_Shift_Del" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'DEL' );
# return "Frame_Shift_Ins" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'INS' );
# return "Nonstop_Mutation" if( $effect eq 'stop_lost' );
# return "Translation_Start_Site" if( $effect =~ /^(initiator_codon_variant|start_lost)$/ );
# return "In_Frame_Ins" if( $effect =~ /^(inframe_insertion|disruptive_inframe_insertion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'INS' ));
# return "In_Frame_Del" if( $effect =~ /^(inframe_deletion|disruptive_inframe_deletion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'DEL' ));
# return "Missense_Mutation" if( $effect =~ /^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)$/ );
# return "Intron" if ( $effect =~ /^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)$/ );
# return "Splice_Region" if( $effect eq 'splice_region_variant' );
# return "Silent" if( $effect =~ /^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$/ );
# return "RNA" if( $effect =~ /^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$/ );
# return "5'UTR" if( $effect =~ /^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$/ );
# return "3'UTR" if( $effect eq '3_prime_UTR_variant' );
# return "IGR" if( $effect =~ /^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)$/ );
# return "5'Flank" if( $effect eq 'upstream_gene_variant' );
# return "3'Flank" if ( $effect eq 'downstream_gene_variant' );

recode_variant_classification <- function(x){
  
  case_when(
    x %in% c(" ", ".", "unknown", "intergenic_variant", "non_coding_transcript_exon_variant") ~ "unknown",
    
    x %in% c("synonymous SNV", "synonymous", "Silent", "synonymous_variant", "stop_retained_variant") ~ "Silent",
    
    x %in% c("3'UTR", "3_prime_UTR_variant", "downstream_gene_variant") ~ "3'UTR",
    
    x %in% c("5'Flank", "upstream_gene_variant", "5_prime_UTR_variant") ~ "5'Flank",
    
    x %in% c("Intron", "intron_variant") ~ "Intron",
    
    # RED
    x %in% c("nonsynonymous SNV", "nonsynonymous", "Missense_Mutation", "NON_SYNONYMOUS_CODING", "missense_variant") ~ "Missense_Mutation",
    
    x %in% c("frameshift deletion", "Frame_Shift_Del") ~ "Frame_Shift_Del",
    x %in% c("nonframeshift deletion", "In_Frame_Del") ~ "In_Frame_Del",
    
    x %in% c("frameshift insertion", "frameshift substitution") ~ "Frame_Shift_Ins",
    x %in% c("In_Frame_Ins") ~ "In_Frame_Ins",
    
    # non-sense
    x %in% c("stoploss SNV", "stoploss", "Nonstop_Mutation", "stop_lost") ~ "Nonstop_Mutation",
    x %in% c("stopgain SNV", "stopgain", "STOP_GAINED", "stop_gained") ~ "Nonsense_Mutation",

    # splicing
    x %in% c("splicing", "splicing*", "Splice_Site", "Splice_Region", "splice_region_variant", "splice_acceptor_variant", "splice_donor_variant") ~ "Splice_Region",
  
    x %in% c("start_lost") ~ "Translation_Start_Site")
}

recode_variant_classification.pcgr <- function(consequence) {
  consequence2 <- case_when(
    consequence == "splice_region_variant&intron_variant" ~ "splice_region_variant",
    consequence == "stop_gained&splice_region_variant" ~ "stop_gained",
    consequence == "missense_variant&splice_region_variant" ~ "splice_region_variant",
    consequence == "splice_region_variant&synonymous_variant" ~ "splice_region_variant",
    TRUE ~ as.character(consequence)
  )
  consequence2 <- recode_variant_classification(consequence2)
  consequence2
}


maf2vcf <- function(maf_fl, 
                    maf2vcf_pl = "",
                    ref_fa = ""
                    
){
  #"module load ensembl_vep_"
  cmd_maf2vcf = glue("{maf2vcf_pl} --input-maf maf_fl --output-dir vcfs --per-tn-vcfs --ref_fasta {ref_fa}")
  
}




# convert mutect1output to maf
# then use maf2vcf to convert into vcf
# useful for ALL CGL output

# ** example mutec2maf --------------
if(FALSE){
  
  library(pacman)
  p_load(dplyr, readr, janitor, magrittr)
  x = "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/mutect/WEX-334187-T__WEX-334187-N_merged.muTect_call_stats.txt"
  chrom = "contig"
  start_pos = "position"
  end_pos = "position"
  
  ref_allele = "ref_allele"
  alt_allele = "alt_allele"
  
  filter = "judgement"
  
  # AF
  t_alt_count = "t_alt_count"
  t_ref_count = "t_ref_count"
  
  sample_name = "tumor_name"
  ref_name = "sample_name"
  
  # mut qual
  mutation_score = "t_lod_fstar"
  
  
}


mutect2maf <- function(x,
                       
                       chrom = "contig",
                       start_pos = "position",
                       end_pos = "position",
                       
                       ref_allele = "ref_allele",
                       alt_allele = "alt_allele",
                       
                       filter = "judgement",
                       
                       # AF
                       t_alt_count = "t_alt_count",
                       t_ref_count = "t_ref_count",
                       
                       sample_name = "tumor_name",
                       ref_name = "sample_name",
                       
                       # mut qual
                       mutation_score = "t_lod_fstar"
                       
                       
                       # TCGA/clinVAR etc
                       
                       
){
  
  if(is.data.frame(x))
    x = as.data.frame(x, stringsAsFactors = FALSE)
  else
    x = read_tsv(x, col_types = cols(.default = col_character())) %>% 
      data.frame(stringsAsFactors = FALSE)
  
  # p_load(GenomicRanges, VariantAnnotation)
  # gr = mutate(x,
  #        end = position) %>% 
  #   dplyr::select(chr = contig, 
  #          start = position,
  #          end,
  #          ref_allele, alt_allele) %>% DataFrame() %>% GRanges()
  # vcf.df = as(elementMetadata(vcf), "DataFrame")
  # anno.df = as(anno, "DataFrame")
  # elementMetadata(vcf) = cbind(vcf.df, anno.df)
  # 
  
  
  #cols_data = colnames(x)
  cols_expected = c(chrom, start_pos, end_pos, 
                    func, ref_allele, alt_allele, 
                    t_alt_count, t_ref_count, 
                    sample_bam, ref_name, sample_name, 
                    mutation_score, aa_change)
  expect_columns(x, cols_expected)
  # columns not included:
  
  # conv imp columns into integer
  x[, t_alt_count] = as.integer(x[, t_alt_count])
  x[, t_ref_count] = as.integer(x[, t_ref_count])
  
  message("> variant type")
  x$.variant_type = get_variant_type(x[, "ref_allele"], x[, "alt_allele"])
  
  message("> seq alleles")
  Tumor_Seq_Allele1 = apply(x, 1, getTumorRef,
                            what = 1,
                            trCount = t_ref_count,
                            rAllele = ref_allele,
                            aAllele = alt_allele, 
                            variant_type = ".variant_type")
  Tumor_Seq_Allele2 = apply(x, 1,  getTumorRef, what = 2,
                            trCount = t_ref_count,
                            rAllele = ref_allele,
                            aAllele = alt_allele, 
                            variant_type = ".variant_type")
  Reference_Allele = get_reference_allele(x[, ref_allele], x$.variant_type)
  
  
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
                    Variant_Type = x$.variant_type,
                    
                    Reference_Allele = Reference_Allele,
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
                    
                    Protein_Change = x[, aa_change],
                    
                    
                    stringsAsFactors = FALSE)
  #categ = apply(x, 1, getCateg, context=context, aAllele=alt_allele, rAllele=ref_allele))
  return(maf)
  
}
mutect_ann_to_maf = tsv2maf






# convert mutect1output to maf
# then use maf2vcf to convert into vcf
# useful for ALL CGL output

# extra --------
# use recode_variant_classification instead
getEffect <- function(mut,
                      func = "func",
                      silent = c(
                        "", NA, "NA", NULL,
                        "unknown", "synonymous SNV",
                        "splicing synonymous SNV"
                      ),
                      nonsilent = c(
                        "nonsynonymous SNV",
                        "stopgain", "stoploss"
                      ),
                      noncoding = ".") {
  if (mut[func] %in% silent) {
    return("silent")
  } else if (mut[func] %in% nonsilent) {
    return("nonsilent")
  } else if (mut[func] %in% noncoding) {
    return("noncoding")
  } else {
    return("null")
  }
}
