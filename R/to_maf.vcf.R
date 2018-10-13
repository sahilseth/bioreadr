
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

#' to maf
#'
#' @param df is a output of parse_somatic_vcf
#'
to_maf.vcf <- function(df){
}