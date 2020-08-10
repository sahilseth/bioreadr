

get_yend <- function(job_no){
  len = length(job_no)
  if(len < 5)
    return(3)
  else
    return(0.5)
}


plot_gantt <- function(x){
  
  #require(ggthemes)
  require(RColorBrewer)
  require(reshape2)
  require(ggplot2)
  require(scales)
  
  x = filter(x, variable %in% c("bgn_time", "end_time")) %>%
    select(jobname, variable, value, everything())
  x$value = norm_offset(x$value)
  x = dcast(x, formula = wd + jobname + job_no + job_id~ variable, value.var = "value") %>% tbl_df()
  
  #   x$bgn_time = as.POSIXct(strptime(x$bgn_time, "%Y-%m-%d %H:%M:%S"))
  #   x$end_time = as.POSIXct(strptime(x$end_time, "%Y-%m-%d %H:%M:%S"))
  
  
  x = group_by(x, jobname) %>%
    mutate(size = as.numeric(get_yend(job_no)))
  
  ##                plot
  p <- ggplot(x) +
    #geom_rect(aes(xmin=bgn_time, xmax=end_time, ymin=job_no, ymax=job_no+size, fill = wd)) +
    geom_rect(aes(xmin=bgn_time, xmax=end_time, ymin=job_no, ymax=job_no+size)) +
    facet_grid(jobname~., scales = "free_y", space = "free") +
    #facet_wrap(~jobname, ncol = 1, shrink = TRUE) +
    #theme_bw() +
    ylab("") + xlab("hours") +
    scale_x_datetime(
      breaks = date_breaks("1 hours"),
      minor_breaks = date_breaks("30 min"),
      labels = date_format("%H:%M")) +
    scale_fill_manual(values=rep(brewer.pal(8, "Set1"), times = 4), guide = FALSE)+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      #axis.title.x	 = element_text(size = 18),
      axis.title.x	 = element_blank(),
      #axis.text.x	 = element_text(size = 18),
      #strip.text.y = element_text(angle = 0, size = 14, face = "bold", color = "gray20"),
      strip.text.y = element_text(angle = 0, color = "gray20"),
      strip.background = element_rect(colour = "gray80", fill = "white"),
      panel.background = element_rect(fill = "gray98"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(size = 0.5, color = "gray85"),
      panel.grid.minor.x = element_line(size = 0.2, color = "gray85")
    )
  
  p
}




get_yend <- function(job_no){
  len = length(job_no)
  if(len < 5)
    return(3)
  else
    return(0.5)
}


get_label_pos <- function(x, y, where = "mid"){
  if(where == "mid"){
    z = (as.numeric(x) + as.numeric(y))/2
  }else if(where == "left"){
    z = as.numeric(y) + 300
  }
  
  as.POSIXct(z, origin = "1970-01-01")
}



#' Title
#'
#' @param x time in character format
#'
norm_offset <- function(x, format = "%Y-%m-%d %H:%M:%S",
                        origin = "1970-01-01"){
  if(is.character(x))
    x = as.POSIXct(strptime(x, format))
  
  offset = min(x) - round_date(min(x), "day")
  
  x = x - offset
  
  ## push back in POSIXct format
  as.POSIXct(x, origin = origin)
}

plot_gantt2 <- function(x, 
                        summarize_jobs = TRUE,
                        pdf = TRUE,
                        colors = c("red", "green", "black"),
                        color_coloumn = "host",
                        fill_coloumn = "host",
                        origin = "1970-01-01", 
                        nudge_y_text = 0.5,
                        cex_text = 3
                        
                        ){
  
  #require(ggthemes)
  require(RColorBrewer)
  require(reshape2)
  require(ggplot2)
  require(scales)
  
  
  x = filter(x, variable %in% c("bgn_time", "end_time")) %>%
    select(jobname, variable, value, everything())
  x$value = norm_offset(x$value)
  
  # since we need to add a arbotrary variable, need to 
  # use paste and formula to add it
  # if both columns are same, dcast becomes erronous
  # thus need to add only unique variables
  # this may throw a error if any of the other vars are repeated.
  my_formula = sprintf("wd + jobname + jobnm + job_no + job_id + %s ~ variable",
                       paste(unique(c(color_coloumn, fill_coloumn)), collapse = "+"))
  x = dcast(x, formula = my_formula,
            value.var = "value") %>% tbl_df()
  
  ## need to transform again, dcast destroys this
  #x$bgn_time = as.POSIXct(strptime(x$bgn_time, "%Y-%m-%d %H:%M:%S"))
  #x$end_time = as.POSIXct(strptime(x$end_time, "%Y-%m-%d %H:%M:%S"))
  x$bgn_time = as.POSIXct(x$bgn_time, origin = origin)
  x$end_time = as.POSIXct(x$end_time, origin = origin)
  
  
  x = group_by(x, jobname) %>%
    mutate(size = as.numeric(get_yend(job_no)))
  
  ## include all the columns we need here
  if(summarize_jobs){
    x = group_by_(x, 'jobnm', 'jobname', fill_coloumn, color_coloumn) %>%
      summarise(size = 0.5,
                bgn_time = min(bgn_time),
                end_time = max(end_time)) %>%
      ungroup() %>%
      arrange(jobname)
    
    x$job_no = nrow(x):1
  }
  
  x_min = min(x$bgn_time)
  x_max = max(x$end_time)
  
  ## get plot aesthetics
  if(pdf){
    line_size = 1
    point_size = 2
  }
  
  #vadjust = 0.6
  ##                plot
  p <- ggplot(x, aes_string(x = 'bgn_time', 
                            xend ='end_time', 
                            y = 'job_no', 
                            yend ='job_no', 
                            color = color_coloumn, 
                            fill = fill_coloumn)) +
    geom_segment(size=line_size, alpha = 0.8) +
    geom_point(aes(x=bgn_time), size=point_size, alpha = 0.8) +
    geom_point(aes(x=end_time), size=point_size, alpha = 0.8) +
    scale_color_manual(values = colors, guide = FALSE) + 
    scale_fill_manual(values = colors, guide = FALSE)
  
  
  ## fix labels
  p <- p + geom_text(aes(x = get_label_pos(bgn_time, end_time), 
                         label = jobnm), 
                     cex = cex_text,  nudge_y = nudge_y_text)
  #p = direct.label(p, "first.qp")
  
  p = p + scale_x_datetime(limits = c(x_min, x_max))
  #fontface = "bold")
  
  p <- p + theme_minimal()
  
  ## remove labels
  p + theme(
    axis.ticks.y = element_blank(),
    axis.text.y	 = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank())
  
  ## 
}



plot_resources <- function(dat){
  #   fmt <- function(){
  #     function(x) format(x,nsmall = 2,scientific = FALSE)
  #   }
  
  mytheme <- theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), 
          axis.title.x = element_blank())
  
  dat2 = tbl_df(dat) %>% 
    filter(variable %in% c('avg_mem', 'max_mem', 'max_swap', 'cpu_time' )) %>%
    mutate(trimmed_value = format(as.numeric(value), digits = 2))
  
  #dat2 = filter(dat2, variable == "cpu")
  #str(dat2)
  
  
  p <- ggplot(dat2, aes(x = jobname, y = as.numeric(trimmed_value))) +
    geom_boxplot() + geom_jitter(col = "grey", alpha = 0.3) + 
    mytheme + ylab("memory in mega bytes; time in hours") +
    facet_wrap(~variable, scales = "free_y")
  
  return(p)
}


plot_time <- function(x, pdf = TRUE, ...){
  
  library(lubridate)
  library(tidyr)
  #x = dat
  x = filter(x, variable %in% c("bgn_time", "end_time"))
  x$value = norm_offset(x$value)
  #x$value = as.POSIXct(strptime(x$value, "%Y-%m-%d %H:%M:%S"))
  x2 = spread(x, variable, value)
  #x2 = dcast(x, wd + jobname + job_id + job_no ~ variable, value.var = "value")
  #x2$interval = interval(x2$bgn_time, x2$end_time)
  x2 = select(x2, jobname, bgn_time, end_time, cores)
  
  x2$cores = as.numeric(x2$cores)
  
  x_min = min(as.POSIXct(x2$bgn_time, origin = origin))
  x_max = max(as.POSIXct(x2$end_time, origin = origin))
  
  ## start from first, and make 1 minute intervals
  time_range <- seq(min(x2$bgn_time), max(x2$end_time), 30) ## total
  cores <- sapply(time_range, function(t){
    rows = which(x2$bgn_time <= t & x2$end_time >= t)
    cpu = sum(x2$cores[rows])
  })
  
  df = data.frame(time = as.POSIXct(time_range, origin = origin), 
                  cores = cores, stringsAsFactors = FALSE)
  
  
  ## plot the result
  p <- ggplot(df) + geom_ribbon(aes(x = time, ymin=0, ymax=cores), alpha = 0.6)
  
  p = p + scale_x_datetime(limits = c(x_min, x_max))
  
  p + theme_minimal() +
    theme(
      axis.ticks.x = element_blank(), 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      
      #axis.ticks.y = element_blank(), 
      #axis.text.y = element_blank(), 
      axis.title.y = element_blank())
  
}

#' Create a plot for gantt and time
#'
#' @param dat 
#'
#' @export
#' 
plot_time_gantt <- function(dat, ...){
  library(grid)
  
  pt = plot_time(dat, ...)
  pg2 = plot_gantt2(dat, ...)
  
  print(pt)
  print(pg2)
  
  grid.newpage()
  grid.draw(rbind(ggplotGrob(pt), ggplotGrob(pg2), size = "first"))
  
  grid.newpage()
  layout = grid.layout(2,1, heights = c(0.2, 0.80))
  pushViewport(viewport(layout = layout))
  print(pt, vp=viewport(layout.pos.row=1, layout.pos.col=1))
  print(pg2, vp=viewport(layout.pos.row=2, layout.pos.col=1))
  
}

#' Title
#'
#' @param x 
#' @param outfile 
#'
#' @export
#'
plot_all <- function(x, outfile){
  library(dplyr)
  library(params)
  library(flowr)
  library(ggplot2)
  library(grid)
  
  set_opts(verbose = 2)
  
  if(missing(outfile))
    outfile = basename(x)
  
  wds = get_wds(x)
  
  fobj = flowr:::read_fobj(wds[1])
  plot_flow(fobj, pdffile = paste0(outfile, "_def.pdf"))
  
  if(!file.exists(paste0(outfile, ".tsv"))){
    
    lst = lapply(wds, function(x){
      flowr:::get_resources_lsf(x)
    })
    
    dat = bind_rows(lst)
    message("writing resources details to file...")
    write_sheet(dat, paste0(outfile, ".tsv"))
    
  }else{
    message("resource details available, reading them...")
    dat = read_sheet(paste0(outfile, ".tsv"))
  }
  
  
  dat = tbl_df(dat)
  message("plotting...")
  pdf(paste0(outfile, ".pdf"))
  
  #pg1 = plot_gantt(dat);print(pg1)
  pr = plot_resources(dat);print(pr)
  
  plot_time_gantt(dat)  
  
  dev.off()
}

