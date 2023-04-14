# libraries
x1 <- c("here","tidyverse","lubridate","Rcpp","RcppRoll","wpp2019") 
# "deSolve","gtools","pracma","lhs",,"ISOweek","tidybayes", "rstudioapi","devtools","wpp2019","epiR"
x2 <- x1 %in% row.names(installed.packages())
if (any(x2 == FALSE)) { install.packages(x1[!x2],repos="http://cran.us.r-project.org") }
# Load all packages (unless already loaded) # as.Date <- zoo::as.Date
lapply(x1,library,character.only=TRUE)
select <- dplyr::select; # row_number <- dplyr::row_number; summarise <- dplyr::summarise

# theme for plots
standard_theme=theme(# panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
  plot.title=element_text(hjust=0.5,size=16),
  axis.text.x=element_text(size=9,angle=90,vjust=1/2),axis.text.y=element_text(size=9),
  axis.title=element_text(size=14), text=element_text(family="Calibri"))


### assign multiple variables ----------------------------------------------------------
# use as: g(a,b,c) %=% c(1,2,3)
'%=%' = function(l, r, ...) UseMethod('%=%')
# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}
# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

### for plots with factors as x-axis, show every nth tick --------------
show_every_nth = function(n) { return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]}) }

### get objects larger than x Mb (memory) --------------
fcn_objs_mem_use <- function(min_size){
  mem_use_df=round(data.frame(unlist(sapply(ls(envir=.GlobalEnv), function(n) object.size(get(n)), simplify = FALSE)))/1e6,1)
  colnames(mem_use_df)[1]<-"size (Mb)"; mem_use_df[,"objs"]=rownames(mem_use_df)
  mem_use_df<-mem_use_df[order(mem_use_df$size,decreasing=T),]; rownames(mem_use_df)<-c()
  mem_use_df[mem_use_df$size>min_size,c(2,1)]
}
