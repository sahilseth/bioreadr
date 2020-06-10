

# Annotating your VCF with VEP
# Adding coverage data to your VCF
# Adding expression data to your VCF
# Creating a phased VCF of proximal variants

# maf>VCF-VEP
# addressed using esembl_vep


# adding read counts -------
# https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/readcounts.html
# Running the vcf-readcount-annotator
# If you have multiple files for SNVs and InDels you will first need to concatenate the two files together:
#   
#   cat snvs_bam_readcount_file indels_bam_readcount_file > bam_readcount_file
# You can now use the combined bam-readcount output file to add readcount information to your VCF.
# 
# vcf-readcount-annotator input_vcf bam_readcount_file DNA|RNA -s sample_name
# vcf-expression-annotator input_vcf expression_file kallisto|stringtie|cufflinks|custom gene|transcript


# phased variants -------
# cases where there are multiple variations in proximity
# https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/proximal_vcf.html
# Update sample names
# The sample names in the tumor.bam, the somatic.vcf, and the germline.vcf need to match. 
# If they don’t you need to edit the sample names in the VCF files to match the tumor BAM file.

# Combine somatic and germline variants using GATK’s CombineVariants
# /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
# -T CombineVariants \
# -R reference.fa \
# --variant germline.vcf \
# --variant somatic.vcf \
# -o combined_somatic_plus_germline.vcf \
# --assumeIdenticalSamples
# Sort combined VCF using Picard
# /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar SortVcf \
# I=combined_somatic_plus_germline.vcf \
# O=combined_somatic_plus_germline.sorted.vcf \
# SEQUENCE_DICTIONARY=reference.dict
# Phase variants using GATKs ReadBackedPhasing
# /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar \
# -T ReadBackedPhasing \
# -R reference.fa \
# -I tumor.bam \
# --variant combined_somatic_plus_germline.sorted.vcf \
# -L combined_somatic_plus_germline.sorted.vcf \
# -o phased.vcf
# bgzip and index the phased VCF
# bgzip -c phased.vcf > phased.vcf.gz
# tabix -p vcf phased.vcf.gz
# The resulting phased.vcf.gz file can be used as the input to the --phased-proximal-variants-vcf option.
# 
# bgzip and index the input VCF
# In order to use the --phased-proximal-variants-vcf option you will also need to bgzip and index your main input VCF.
# 
# bgzip -c input.vcf > input.vcf.gz
# tabix -p vcf input.vcf.gz


# pvacseq --------

#' @name pvacseq
#'
#' @param input_vcf input_vcf
#' @param samplename something
#' @param hla_alleles  something
#' @param mhc_tools something
#' @param pvacseq_opts something
#' @param odir something
#' @param pvacseq_setup something
#' 
#' @details Still only works on einstein
#'
#' @export
pvacseq <- function(
  input_vcf,
  hla_fls,
  #hla_alleles = "HLA-A*02:01,HLA-B*35:01,DRB1*11:01",
  samplename,
  seq_sample_id_t,
  
  odir = "pvacseq_full",
  
  # pvacseq_mhc_tools = "MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign",
  pvacseq_mhc_tools = opts_flow$get("pvacseq.mhc_tools"),
  # pvacseq_mhc_tools = "MHCnuggetsI MHCnuggetsII",
  # ;source activate pvactools;export KERAS_BACKEND=tensorflow
  # pvacseq_exe = "module load pvactools_/1.5.8",
  pvacseq_setup = opts_flow$get("pvacseq.setup"),
  pvacseq_opts = opts_flow$get("pvacseq.opts"),
  
  iedb_dir = opts_flow$get("pvacseq.iedb_dir")
  #"/home/sseth/apps/iedb"
  
){
  check_args()
  
  # pvacseq run ${ann_vcf} "1004-t" $(cat WEX-1004-N_nochr_hla_umap_hla_pvac.txt) ${mhc_tools} pvacseq -e 8,9,10 --iedb-install-directory $iedb_dir
  # skipping iedb for now, causes major issues
  # https://unix.stackexchange.com/questions/37135/concatenate-two-files-without-adding-a-newline
  cmd_pvacseq = glue("{pvacseq_setup}; ",
                     "rm -rf {odir}; ",
                     "pvacseq run {input_vcf} {seq_sample_id_t} $(paste -d',' {hla_fls}) {pvacseq_mhc_tools} {odir} {pvacseq_opts} --iedb-install-directory {iedb_dir}")
  
  flowmat = to_flowmat(list(pvacseq = cmd_pvacseq), samplename)
  
  outfiles = list(epitopes_all = glue("pvacseq/MHC_Class_I/{seq_sample_id_t}.all_epitopes.tsv"),
                  epitopes_filt = glue("pvacseq/MHC_Class_I/{seq_sample_id_t}.filtered.tsv"),
                  epitopes_rnk = glue("pvacseq/MHC_Class_I/{seq_sample_id_t}.filtered.condensed.ranked.tsv"))
  
  list(flowmat = flowmat, outfiles = outfiles)
}

pvac_prep_vcf <- function(input_vcf = input_vcf, 
                          samplename = samplename,
                          
                          vcf_ext = "vcf",
                          seq_sample_id_wex_t,
                          seq_sample_id_wex_n,
                          seq_sample_id_rna_t,
                          seq_sample_id_wex_t_vcf,
                          
                          bam_rna_t = NA,
                          exp_rna_t = NA,
                          exp_rna_type = c("rsem", "kallisto_gene"),
                          
                          ref_fasta = opts_flow$get("ref_fasta"),
                          bam_readcount_helper_py = opts_flow$get("bam_readcount_helper_py"),
                          
                          odir = "pvac_prep_vcf"
                          
                          
){
  
  exp_rna_type = match.arg(exp_rna_type)
  check_args()
  
  ann_vcf = gsub(".vcf$|.vcf.gz$", "_vep.vcf", input_vcf)
  ann_vcf = file.path(odir, ann_vcf)
  out_vep = ensembl_vep(input_vcf = input_vcf,
                        samplename = samplename,
                        ann_vcf = ann_vcf)
  curr_vcf = ann_vcf
  cmds = out_vep$flowmat$cmd
  cmds = glue("mkdir {odir};{cmds}")
  
  if(!is.na(bam_rna_t)){
    # output_vcf="IPCT-S4008-MOON0051-Cap2043-4-ID01_190415-A00422-0045-AHK77HDSXX-4-ATCACG--S4008-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG_vep_decomp.vcf"
    decomp_vcf = gsub(".vcf$", ".decomp.vcf", curr_vcf)
    cmd_decomp = glue("vt decompose -s {ann_vcf} -o {decomp_vcf}")
    curr_vcf = decomp_vcf
    
    # run bam read count
    # # <decomposed_vcf> <sample_name> <reference_fasta> <bam_file> <output_dir>
    cmd_readcount = glue("mkdir bam_readcount_rna;{bam_readcount_helper_py} {curr_vcf} {seq_sample_id_wex_t_vcf} {ref_fasta} {bam_rna_t} bam_readcount_rna")
    
    # annotate vcf
    bam_readcount_rna_snv = glue("bam_readcount_rna/{seq_sample_id_wex_t_vcf}_bam_readcount_snv.tsv")
    # echo $bam_readcount_rna_snv
    cmd_vcf_ann_rna_vaf = glue("vcf-readcount-annotator {curr_vcf} {bam_readcount_rna_snv} RNA -s {seq_sample_id_wex_t_vcf}")
    readcount_vcf = gsub(".vcf$", ".readcount.vcf", curr_vcf)
    curr_vcf = readcount_vcf
    
    cmds = c(cmds, 
             cmd_decomp, 
             cmd_readcount, 
             cmd_vcf_ann_rna_vaf)
  }
  
  # https://vatools.readthedocs.io/en/latest/vcf_expression_annotator.html
  
  # The column header in the expression_file for the
  # column containing gene names/transcript ids. Required
  # when using the `custom` format.
  # for Kallisto it uses symbols, for all others it uses gene_ids!
  # if args.format == 'kallisto':
  #     if key == 'SYMBOL' and value != '':
  #         genes.add(value)
  # else:
  #     if key == 'Gene' and value != '':
  #         genes.add(value)
  ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|DownstreamProtein|ProteinLengthChange|WildtypeProtein">
  # "CSQ=A|splice_region_variant&intron_variant|LOW|NPHP4|ENSG00000131697|Transcript|ENST00000378156|protein_coding||27/29||||||||||-1||HGNC|19104|||MNDWHRIFTQNVLVPPHPQRARQPWKESTAFQCVLKWLDGPVIRQGVLEVLSEVECHLRVSFFDVTYRHFFGRTWKTTVKPTKRPPSRIVFNEPLYFHTSLNHPHIVAVVEVVAEGKKRDGSLQTLSCGFGILRIFSNQPDSPISASQDKRLRLYHGTPRALLHPLLQDPAEQNRHMTLIENCSLQYTLKPHPALEPAFHLLPENLLVSGLQQIPGLLPAHGESGDALRKPRLQKPITGHLDDLFFTLYPSLEKFEEELLELHVQDHFQEGCGPLDGGALEILERRLRVGVHNGLGFVQRPQVVVLVPEMDVALTRSASFSRKVVSSSKTSSGSQALVLRSRLRLPEMVGHPAFAVIFQLEYVFSSPAGVDGNAASVTSLSNLACMHMVRWAVWNPLLEADSGRVTLPLQGGIQPNPSHCLVYKVPSASMSSEEVKQVESGTLRFQFSLGSEEHLDAPTEPVSGPKVERRPSRKPPTSPSSPPAPVPRVLAAPQNSPVGPGLSISQLAASPRSPTQHCLARPTSQLPHGSQASPAQAQEFPLEAGISHLEADLSQTSLVLETSIAEQLQELPFTPLHAPIVVGTQTRSSAGQPSRASMVLLQSSGFPEILDANKQPAEAVSATEPVTFNPQKEESDCLQSNEMVLQFLAFSRVAQDCRGTSWPKTVYFTFQFYRFPPATTPRLQLVQLDEAGQPSSGALTHILVPVSRDGTFDAGSPGFQLRYMVGPGFLKPGERRCFARYLAVQTLQIDVWDGDSLLLIGSAAVQMKHLLRQGRPAVQASHELEVVATEYEQDNMVVSGDMLGFGRVKPIGVHSVVKGRLHLTLANVGHPCEQKVRGCSTLPPSRSRVISNDGASRFSGGSLLTTGSSRRKHVVQAQKLADVDSELAAMLLTHARQGKGPQDVSRESDATRRRKLERMRSVRLQEAGGDLGRRGTSVLAQQSVRTQHLRDLQVIAAYRERTKAESIASLLSLAITTEHTLHATLGVAEFFEFVLKNPHNTQHTVTVEIDNPELSVIVDSQEWRDFKGAAGLHTPVEEDMFHLRGSLAPQLYLRPHETAHVPFKFQSFSAGQLAMVQASPGLSNEKGMDAVSPWKSSAVPTKHAKVLFRASGGKPIAVLCLTVELQPHVVDQVFRFYHPELSFLKKAIRLPPWHTFPGAPVGMLGEDPPVHVRCSDPNVICETQNVGPGEPRDIFLKVASGPSPEIKDFFVIIYSDRWLATPTQTWQVYLHSLQRVDVSCVAGQLTRLSLVLRGTQTVRKVRAFTSHPQELKTDPKGVFVLPPRGVQDLHVGVRPLRAGSRFVHLNLVDVDCHQLVASWLVCLCCRQPLISKAFEIMLAAGEGKGVNKRITYTNPYPSRRTFHLHSDHPELLRFREDSFQVGGGETYTIGLQFAPSQRVGEEEILIYINDHEDKNEEAFCVKVIYQ"
  # NPHP4|ENSG00000131697|Transcript|ENST00000378156
  # exp:
  # ENST00000489180.6|ENSG00000131697.17|OTTHUMG00000000701.6|OTTHUMT00000001716.4|NPHP4-209|NPHP4|5548|UTR5:1-268|CDS:269-3004|UTR3:3005-5548|     5548    5349    32.1671 0.523718
  # ENST00000378156.8|ENSG00000131697.17|OTTHUMG00000000701.6|OTTHUMT00000001715.2|NPHP4-201|NPHP4|4994|UTR5:1-266|CDS:267-4547|UTR3:4548-4994|     4994    4795    144.317 2.62111
  # ENST00000378169.7|ENSG00000131697.17|OTTHUMG00000000701.6|OTTHUMT00000367612.1|NPHP4-203|NPHP4|4213|UTR5:1-242|CDS:243-1031|UTR3:1032-4213|     4213    4014    47.0472 1.02074
  # ENST00000622020.4|ENSG00000131697.17|OTTHUMG00000000701.6|-|NPHP4-211|NPHP4|3004|UTR5:1-268|CDS:269-3004|       3004    2805    20.5988 0.639539
  # ENST00000466897.1|ENSG00000131697.17|OTTHUMG00000000701.6|OTTHUMT00000001719.3|NPHP4-205|NPHP4|1192|CDS:1-415|UTR3:416-1192|    1192    993     17.8703 1.56725
  
  # example data: 
  # https://github.com/griffithlab/VAtools/blob/e500e92e801a06ed31b7459466ffc1a21a341c6a/tests/test_data/kallisto.genes
  # gene_name	gene	abundance	counts	length
  # A1BG	A1BG	9.5072780253347	253.812537262238	1113.77371499879
  # A1CF	A1CF	0.170005984522626	33.3184303003113	8176.36764900294
  # A2M	A2M	4.49761004558132	110	1020.35424350425
  if(!is.na(exp_rna_t)){
    if(exp_rna_type == "rsem"){
      cmd_vcf_ann_rna_exp = glue("vcf-expression-annotator -e 'TPM' -i 'gene_id' -s {seq_sample_id_wex_t_vcf} {curr_vcf} {exp_rna_t} custom gene")
    }else if(exp_rna_type == "kallisto"){
      cmd_vcf_ann_rna_exp = glue("vcf-expression-annotator -s {seq_sample_id_wex_t_vcf} {curr_vcf} {exp_rna_t} kallisto transcript")
    }else if(exp_rna_type == "kallisto_gene"){
      cmd_vcf_ann_rna_exp = glue("vcf-expression-annotator -e 'tpm' -i 'gene_id' -s {seq_sample_id_wex_t_vcf} --ignore-ensembl-id-version {curr_vcf} {exp_rna_t} custom gene")
    }
    
    gx_vcf = gsub(".vcf$", ".gx.vcf", curr_vcf)
    curr_vcf = gx_vcf
    cmds = c(cmds, cmd_vcf_ann_rna_exp)
  }
  
  
  flowmat = to_flowmat(list(pvac_prep_vcf = cmds), samplename)
  
  list(flowmat = flowmat, outfiles = list(final_vcf = curr_vcf))
  
  
}

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
