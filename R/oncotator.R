
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
#' @param x tsv file, Hugo_Symbol, Chromosome, Start_Position, End_Position,Reference_Allele, Tumor_Seq_Allele2,
#' @param oncotator_exe 
#' @param db_path 
#' @param outfile 
oncotator <- function(x, outfile,
                      oncotator_exe = "/home/sseth/apps/python/2.7/bin/oncotator",
                      db_path = "~/rsrch2_home/ref/human/b37/annotations/oncotator/oncotator_v1_ds_April052016",
                      execute = F){
  
  library(glue)
  
  # a few lines to confirm
  # this is a valid maf format
  if(is.data.frame(x)){
    maflite = tempfile()
    write_tsv(x, maflite)
  }else{
    maflite = tools::file_path_as_absolute(x)
  }
  #outfile = tools::file_path_as_absolute(outfile)
  
  cmd_onco = glue("{oncotator_exe} -v --db-dir {db_path} {maflite} {outfile} hg19")
  print(cmd_onco)
  
  if(execute){
    system(cmd_onco)
    
    df = readr::read_tsv(outfile, skip = 3) %>% clean_names()
    df
  }else{
    cmd_onco
  }
}



# using API -----
#   chromsome    start      end     ref alt Tumor_Sample_Barcode
# 1      chr4 55589774 55589774       A   G               fake_1
# 2      chr4 55599321 55599321       A   T               fake_2
# 3      chr4 55599332 55599332       G   T               fake_3
# 4      chr4 55599320 55599320       G   T               fake_4
# 5     chr15 41961117 41961123 TGGCTAA   -               fake_4
# 6      chr4 55599320 55599320       G   T               fake_5



# END

