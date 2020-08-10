


ensembl_vep <- function(input_vcf,
                        samplename,
                        ann_vcf,
                        
                        vep_exe = opts_flow$get("vep_exe"),
                        vep_plugins_dir = opts_flow$get("vep_plugins_dir"),
                        vep_dir = opts_flow$get("vep_dir"),
                        vep_assembly_version = opts_flow$get("vep_assembly_version"),
                        vep_ref_fa = opts_flow$get("vep_ref_fa")
                        
                        ){
  
  # vep
  cmd_vep = glue("{vep_exe} --input_file {input_vcf} --format vcf --output_file {ann_vcf} ",
                 "--vcf --symbol --terms SO --plugin Downstream --plugin Wildtype --dir_plugins {vep_plugins_dir} --assembly {vep_assembly_version} ",
                 "--fasta {vep_ref_fa} --dir_cache {vep_dir} --offline --pick --force_overwrite")
  
  # #| $vep_filter --format vcf --force_overwrite --filter "(FILTER is PASS)" --output_file $annotated_dir/$sampleID"_annotated_filterd.vcf"
  
  # need the following opts
  # --format vcf
  # --vcf
  # --symbol
  # --plugin Downstream
  # --plugin Wildtype
  # --terms SO
  
  flowmat = to_flowmat(list(vep = cmd_vep), samplename)
  
  list(flowmat = flowmat, outfiles = list(ann_vcf = ann_vcf))
  
}
