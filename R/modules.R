
pkg_topic <- function(package, topic, file = NULL) {
  # Find "file" name given topic name/alias
  if (is.null(file)) {
    topics <- pkg_topics_index(package)
    topic_page <- subset(topics, alias == topic, select = file)$file

    if(length(topic_page) < 1)
      topic_page <- subset(topics, file == topic, select = file)$file

    stopifnot(length(topic_page) >= 1)
    file <- topic_page[1]
  }

  rdb_path <- file.path(system.file("help", package = package), package)

  tools:::fetchRdDB(rdb_path, file)
}

pkg_topics_index <- function(package) {
  help_path <- system.file("help", package = package)

  file_path <- file.path(help_path, "AnIndex")
  if (length(readLines(file_path, n = 1)) < 1) {
    return(NULL)
  }

  topics <- read.table(file_path, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "", quote = "", header = FALSE)

  names(topics) <- c("alias", "file")
  topics[complete.cases(topics), ]
}

parse_rd <- function (x, ...){
  tags <- vapply(x, tag, FUN.VALUE = character(1))
  get_tags <- function(tag) x[tags == tag]
  get_tag <- function(tag) {
    if (tag %in% tags) {
      x[[which(tags == tag)]]
    }
  }
  line_breaks <- tags == "TEXT"
  x <- x[!line_breaks]
  tags <- tags[!line_breaks]
  out <- list()
  out$name <- get_tag("name")
  out$title <- get_tag("title")
  out$aliases <- vapply(get_tags("alias"), as.character, character(1),
                        ...)
  out$keywords <- vapply(get_tags("keyword"), as.character, character(1),
                         ...)
  out$usage <- get_tag("usage")
  out$arguments <- get_tag("arguments")
  if (length(out$arguments)) {
    out$has_args <- TRUE
  }
  out$author <- get_tag("author")
  out$seealso <- get_tag("seealso")
  out$examples <- get_tag("examples")
  sections <- x[!(tags %in% c("name", "title", "alias", "keyword",
                              "usage", "author", "seealso", "arguments", "examples"))]
  out
}

#' Show a list of all available modules in ultraseq
#'
#' @param pkg
#' @param show
#'
#' @export
#'
#' @examples
#' modules()
modules <- function(pkg = "ultraseq", show = TRUE){
  
  funcs = ls(paste0('package:', pkg))

  # for each function check attribute type, if module, return name.
  f = funcs[1]
  
  tmp <- lapply(funcs, function(f){
    func = get(f)
    type = attr(func, 'type', exact = TRUE)
    type = ifelse(is.null(type), 'NA', type)

    #debug(parse_rd)
    #rd =  parse_rd(rd)
    ttl = NA

    if(type == "module")
      ttl = pkg_topic(pkg, f)[[1]][[1]]

    data.frame(name = f, type = type, desc = ttl,
               stringsAsFactors = FALSE)
  })
  
  tb = do.call(rbind, tmp)
  tb = subset(tb, type == "module")



  if(show)
    kable(tb)

}
