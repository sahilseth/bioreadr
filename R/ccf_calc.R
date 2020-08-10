
if(FALSE){
  
  
}

calc_u.ccf <- function(ccf, m){
  ccf*m
}

calc_ccf.u <- function(u, m){
  return (u/m)
}

calc_m.u <- function(u){
  if(u >=1)
    m = round(abs(u), 0)
  if(u < 1)
    m = 1
  m
}


calc_u  <- function(f, p, cn){
  (f * 1/p) * (p*cn + (1-p)*2)
  
}
