

vt_decompose <- function(infl, outfl = tempfile(fileext = ".vcf")){
  vt_exe = "~/.conda/envs/pcgr_py36/bin/vt"
  cmd = glue("{vt_exe} decompose -s {infl} > {outfl}");cmd
  system(cmd)
  outfl
}