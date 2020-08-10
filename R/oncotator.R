
# oncodb=~/rsrch2_home/ref/human/b37/annotations/oncotator/oncotator_v1_ds_April052016
# maf_fl=~/.rsrch1/iacs/iacs_dep/sseth/projects2/ss_drug_devel/shp2/wrangl_01_shp2_ccle/df_ccle_shp2_muts_input.tsv
# 
#
# oncotator -v --db-dir $oncodb $maf_fl exampleOutput.tsv hg19

#' creating a maf file
#' Mandatory fields: 
#' Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, 
#' Tumor_Seq_Allele2, Variant_Classification, Variant_Type and Tumor_Sample_Barcode.
#'
#' @param x tsv file, Hugo_Symbol, Chromosome, Start_Position, End_Position,Reference_Allele, Tumor_Seq_Allele2, OR a VCF file
#' @param oncotator_exe 
#' @param db_dir 
#' @param outfile 
oncotator <- function(x, 
                      samplename,
                      outfile,
                      input_format = "maflite",
                      output_format = "maflite",
                      oncotator_exe = "~/apps/python/2.7/bin/oncotator",
                      db_dir = "~/rsrch2_home/ref/human/b37/annotations/oncotator/oncotator_v1_ds_April052016",
                      execute = F, 
                      oncotator_opts = "-c $tx"){
  
  library(glue)
  
  if(input_format == "maflite"){
    # a few lines to confirm
    # this is a valid maf format
    if(is.data.frame(x)){
      input_fl = tempfile()
      write_tsv(x, maflite)
    }
  }else{
    input_fl = sapply(x, tools::file_path_as_absolute)
    outfile = sapply(outfile, tools::file_path_as_absolute)
  }
  
  #cmd_onco = glue("{oncotator_exe} -v --db-dir {db_path} --output_format {output_format} {input_fl} {outfile} hg19")
  cmd_onco = glue("{oncotator_exe} -v --db-dir {db_dir} --output_format {output_format} {oncotator_opts} -i {input_format} {input_fl} {outfile} hg19")
  
  #print(cmd_onco)
  
  if(execute){
    system(cmd_onco)
    
    if(output_format == "maflite"){
      df = readr::read_tsv(outfile, skip = 3) %>% clean_names()
      return(df)
    }
  }else{
    flowmat = to_flowmat(list(cmd_onco = cmd_onco), samplename = samplename)
    return(flowmat)
  }
  
}

# https://gatkforums.broadinstitute.org/gatk/discussion/4154/howto-install-and-run-oncotator-for-the-first-time
# oncotator=~/apps/conda/2.7/envs/oncotator/bin/oncotator
# db_dir=~/ref/human/b37/annotations/broad-somatic-b37/oncotator/oncotator_v1_ds_April052016
# tx=~/ref/human/b37/annotations/broad-somatic-b37/oncotator/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt
# input_fl=sarco5-T_1208XX_ST1374_073_H09WGADXX___sarco5-N_1208XX_ST1374_073_H09WGADXX_filt2.vcf
# outfile==sarco5-T_1208XX_ST1374_073_H09WGADXX___sarco5-N_1208XX_ST1374_073_H09WGADXX_filt2_onco.vcf
# $oncotator -v --db-dir ${db_dir} --output_format VCF -c $tx -i VCF ${input_fl} ${outfile} hg19 


# using API -----
#   chromsome    start      end     ref alt Tumor_Sample_Barcode
# 1      chr4 55589774 55589774       A   G               fake_1
# 2      chr4 55599321 55599321       A   T               fake_2
# 3      chr4 55599332 55599332       G   T               fake_3
# 4      chr4 55599320 55599320       G   T               fake_4
# 5     chr15 41961117 41961123 TGGCTAA   -               fake_4
# 6      chr4 55599320 55599320       G   T               fake_5



# END

