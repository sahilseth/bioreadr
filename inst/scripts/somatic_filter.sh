function somatic_filter() {
  awk -v MINQUAL="$1" -v SSC_THRES="$2" -v ONLY_SOMATIC="$3" 'BEGIN {NORMAL=10; TUMOR=11; GL_IDX=0;}
  {
  if ($0~"^#") { print ; next; }
  if (! GL_IDX) {
  split($9,fmt,":")
  for (i=1;i<=length(fmt);++i) { if (fmt[i]=="GL") GL_IDX=i }
  }
  split($NORMAL,N,":");
  split(N[GL_IDX],NGL,",");
  split($TUMOR,T,":");
  split(T[GL_IDX],TGL,",");
  LOD_NORM=NGL[1]-NGL[2];
  LOD_TUMOR_HET=TGL[2]-TGL[1];
  LOD_TUMOR_HOM=TGL[3]-TGL[1];
  if (LOD_TUMOR_HET > LOD_TUMOR_HOM) { LOD_TUMOR=LOD_TUMOR_HET }
  else { LOD_TUMOR=LOD_TUMOR_HOM }
  DQUAL=LOD_TUMOR+LOD_NORM;
  if (DQUAL>=SSC_THRES && $NORMAL~"^0/0") {
  $7="PASS"
  $8="SSC="DQUAL";"$8
  print
  }
  else if (!ONLY_SOMATIC && $6>=MINQUAL && $10~"^0/0" && ! match($11,"^0/0")) {
  $8="SSC="DQUAL";"$8
  print
  }
  }' OFS="\t"
}
