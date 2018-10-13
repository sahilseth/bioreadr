
## ----doit1---------------------------------------------------------------
library(XML)
top <- newXMLNode("top")
tvp <- newXMLNode("TVP", parent = top)
time <- newXMLNode("time", "2012-01-01", parent = tvp)
value <- newXMLNode("value", "123", parent = tvp)
comment <- newXMLNode("comment",parent = tvp)
qualifer <-newXMLNode("qualifier", attrs = c(y = 'abc'), parent = comment)
commentText <-newXMLNode("info", attrs = c(y = 'something'), parent = comment)
tvp <- newXMLNode("TVP", parent = top)
time <- newXMLNode("time", "2012-01-02", parent = tvp)
value <- newXMLNode("value", "456", parent = tvp)
tvp <- newXMLNode("TVP", parent = top)
time <- newXMLNode("time", "2012-01-03", parent = tvp)
value <- newXMLNode("value", "789", parent = tvp)
comment <- newXMLNode("comment",parent = tvp)
newXMLNode("qualifier", attrs = c(y = 'efg'), parent = comment)
top

### My answer
require(plyr) ### provides rbind.fill
getDataframe <- function(xml){
 out2 <- xmlSApply(xml,function(x){
 out <- xmlSApply(x, function(y){
  if(length(xmlChildren(y)) >= 1){xmlSApply(y,xmlAttrs)
   }else{xmlValue(y)}})
  as.data.frame(t(unlist(out))) ## rbind.fill likes dataframes
  })
  return(do.call(rbind.fill,out2))
}
getDataframe(top)

## http://stackoverflow.com/questions/18538772/parsing-irregular-xml-in-r/18539869#18539869
xp <- xpathApply(top, "/top/TVP", xpathSApply, ".//*[not(*)]", function(x){
  print(x)
  len=nzchar(xmlValue(x));len
  ret=ifelse(len, xmlValue(x), xmlAttrs(x))
  return(setNames(ret,xmlName(x)))
                 })
library(plyr)
do.call(rbind.fill.matrix, lapply(xp, t))




## ----skipit,fig.width=7, fig.height=6------------------------------------
plot(cars)


