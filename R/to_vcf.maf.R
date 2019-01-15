


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
    x %in% c(" ", ".", "unknown") ~ "unknown",
    
    x %in% c("synonymous SNV", "synonymous", "Silent") ~ "Silent",
    
    x %in% c("3'UTR") ~ "3'UTR",
    
    x %in% c("5'Flank") ~ "5'Flank",
    
    x %in% "Intron" ~ "Intron",
    
    # RED
    x %in% c("nonsynonymous SNV", "nonsynonymous", "Missense_Mutation", "NON_SYNONYMOUS_CODING") ~ "Missense_Mutation",
    
    x %in% c("frameshift deletion", "Frame_Shift_Del") ~ "Frame_Shift_Del",
    
    x %in% c("frameshift insertion", "synonymous", "Silent") ~ "Silent",
    
    x %in% c("nonframeshift deletion", "In_Frame_Del") ~ "In_Frame_Del",
    x %in% c("In_Frame_Ins") ~ "In_Frame_Ins",
    
    x %in% c("stoploss SNV", "stoploss", "Nonstop_Mutation") ~ "Nonstop_Mutation",
    x %in% c("stopgain SNV", "stopgain", "STOP_GAINED") ~ "Nonsense_Mutation",
    x %in% c("In_Frame_Ins") ~ "In_Frame_Ins",
    
    x %in% c("splicing", "splicing*", "Splice_Site", "Splice_Region") ~ "Splice_Region")
  
  
}