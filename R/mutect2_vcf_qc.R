
# more here: http://einstein:2500/notebooks/projects2/ss_sarco/dnaseq/02b_m2_vs_ir.ipynb

gg2density <- function(df, x, y){
  p <- ggplot(df, aes_string(x, y)) +
    stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE) +
    scale_fill_viridis_c() +
    coord_cartesian(expand = FALSE) +
    geom_point(shape = '.', col = 'white', alpha = 0.05, size = 0.1)
  p
}







# END


