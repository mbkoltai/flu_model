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

data(pop)

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

### get objects larger than x Mb (memory) --------------

fcn_shortercntr_names <- function(df) { 
  df_out <- df %>% 
  mutate(country_altern_name=country,
         country_altern_name=case_when(
           country %in% "Bolivia (Plurinational State of)" ~ "Bolivia",
           country %in% "Venezuela (Bolivarian Republic of)" ~ "Venezuela",
           country %in% "United States of America" ~ "USA",
           country %in% "Russian Federation" ~ "Russia",
           country %in% "United Kingdom, England" ~ "England",
           country %in% "Iran (Islamic Republic of)" ~ "Iran",
           # country %in% "Côte d'Ivoire" ~ "Cote d'Iv",
           country %in% "Democratic Republic of the Congo" ~ "DRC",
           country %in% "Central African Republic" ~ "CAR",
           country %in% "Lao People's Democratic Republic" ~ "Laos",
           grepl("Tanzania",country) ~ "Tanzania",
           grepl("Korea",country) ~ "S Korea",
           TRUE ~ country))
  return(df_out)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# function to identify flu seasons

fcn_identify_seasons <- function(df_input,sel_variable,source_varname="ORIGIN_SOURCE",
                                 up_thresh=0.9,low_thresh=0.6,length_lim=9,print_flag=F){
  
  df_inds <- df_input %>% ungroup() %>%  select(country,STRAIN,!!source_varname) %>% unique()
  list_cntr_epid_incl<-list()
  
  # sel_variable 
  if (sel_variable %in% "positivity") {
    df_input <- df_input %>% rename(count=value) %>% mutate(value=count/SPEC_PROCESSED_NB)
  }
  
  for (k_row in 1:nrow(df_inds)) {
    sel_source_varname <- as.character(df_inds[k_row,source_varname])
    
    input_data <- df_input %>% 
      filter(country %in% df_inds$country[k_row] & 
                                        STRAIN %in% df_inds$STRAIN[k_row] & 
                                        !!sym(source_varname) %in% sel_source_varname) %>%
      mutate(flu_peak=quantile(value,probs=up_thresh),over_peak=F,
             over_peak=flu_peak<value,
             flu_included=quantile(value, probs=low_thresh),
             over_inclusion=F, over_inclusion=flu_included<value,
             seq_log=rep(rle(over_inclusion)$length>=length_lim,
                         times=rle(over_inclusion)$length))
    if (print_flag){
          print(paste0(paste0(dim(input_data),collapse = ", "),", ",paste0(df_inds[k_row,],collapse =", ")))
      }
    
    tmp <- rle(input_data$over_inclusion)
    
    # add the sequence number to each
    seq_to_add <- c(); start_seq <- 1
    for(i in 1:length(tmp$lengths)){
      length_run <- tmp$lengths[i]
      tester <- input_data[sum(tmp$lengths[1:i]),"over_inclusion"]
      if(tester == F){ 
        seq_to_add <- c(seq_to_add, rep(0,length_run))
      } else {
        seq_to_add <- c(seq_to_add, rep(start_seq,length_run))
        start_seq <- start_seq + 1
      }
    }
    
    input_data$seq <- seq_to_add; input_data$epidem_inclusion <- 0
    for(j in 1:start_seq){
      if (sum(input_data$over_peak[input_data$seq==j])>0) {
        input_data$epidem_inclusion[input_data$seq==j]=1  }
    }
    
    list_cntr_epid_incl[[k_row]] <- input_data %>% rename(positivity=value,value=count)
  }
  
  return(bind_rows(list_cntr_epid_incl))
}
