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
      mutate(flu_peak=quantile(value,probs=up_thresh), over_peak=F, over_peak=flu_peak<value,
             flu_included=quantile(value, probs=low_thresh), 
                    over_inclusion=F, over_inclusion=flu_included<value) 
    # seq_log=rep(rle(over_inclusion)$length>=length_lim,times=rle(over_inclusion)$length)
    if (print_flag){
          print(paste0( paste0(dim(input_data),collapse = ", "),", ",
                        paste0(df_inds[k_row,],collapse =", ")) )
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
      # any data points over the peak?
      peak_detect <- sum(input_data$over_peak[input_data$seq==j])>0
      # is the selected period at least `length_lim` weeks long?
      length_detect <- sum(input_data$seq==j)>=length_lim
      if (peak_detect&length_detect) {
        input_data$epidem_inclusion[input_data$seq==j]=1  }
    }
    
    input_data <- input_data %>%
      mutate(epid_index = cumsum(epidem_inclusion == 1 & lag(epidem_inclusion, default = 0) == 0) ,
             epid_index=ifelse(epidem_inclusion == 0,0,epid_index))
    
    list_cntr_epid_incl[[k_row]] <- input_data %>% rename(positivity=value,value=count)
  }
  
  return(bind_rows(list_cntr_epid_incl))
}


### ### ### ### ### ### 
fcn_find_bloc_lims <- function(df_cont_data,log_flag=T) {
  df_cont_data %>% group_by(country) %>% 
  mutate(min_val=ifelse(log_flag,1,min(value)),max_val=max(value)) %>%
  group_by(country,metasource,STRAIN) %>% mutate(n_data_points=cumsum(epidem_inclusion==0)) %>%
  filter(epidem_inclusion==1) %>%
  group_by(country,metasource,STRAIN,n_data_points) %>%
  summarize(start_date=min(ISO_WEEKSTARTDATE),end_date=max(ISO_WEEKSTARTDATE),
            min_val=unique(min_val),max_val=max(max_val)) %>%
    mutate(country=ifelse(grepl("King",country),"UK",country)) %>% 
    arrange(country, metasource, STRAIN) %>% mutate(index=row_number())
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# create reduced contact matrix ----------------------------------------------------------
fun_create_red_C_m <- function(C_m_full,model_agegroups,orig_age_groups_duration,orig_age_groups_sizes){
  n_age <- nrow(model_age_groups)
  C_m=matrix(0,nrow=n_age, ncol=n_age)
  rownames(C_m)=model_agegroups$agegroup_name; colnames(C_m)=model_agegroups$agegroup_name
  for (i_row in 1:n_age){
    for (j_col in 1:n_age){
      # we are merging or splitting age groups, there are 3 possibilities for an age group in the NEW matrix:
      # 1) same 2) smaller than the group in original CM 3) larger than group in the original CM
      # (it is an *average* of contacts per person, and we have no resolution within age bands)
      #
      # if the 'i' group (C[i,j]) is the *same* as original or *smaller* contact rate UNCHANGED
      if (model_agegroups$wpp_agegroup_low[i_row]==model_agegroups$wpp_agegroup_high[i_row]) {
        # if contact (j) group same or smaller as original
        if (model_agegroups$wpp_agegroup_low[j_col]==model_agegroups$wpp_agegroup_high[j_col]) {
          # proportionality by 'time window' size (age x to x+n)
          f_dur=model_agegroups$duration[j_col]/orig_age_groups_duration[model_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=(C_m_full[model_agegroups$wpp_agegroup_low[i_row],
                                     model_agegroups$wpp_agegroup_low[j_col]])*f_dur
          # print("'i' group (C[i,j]) is the same/smaller as original & 'j' group same/smaller as original"); 
        } else { 
          # contact numbers BY contact groups (j) we need to SUM
          group_span=model_agegroups$wpp_agegroup_low[j_col]:model_agegroups$wpp_agegroup_high[j_col]
          agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
          C_m[i_row,j_col]=sum(C_m_full[i_row,group_span]) # agegroup_weights*
          # print("'i' group (C[i,j]) is the *same* as original or *smaller* & 
          #       'j' is larger than original group"); print(c(i_row,j_col))
        } # end of 'i' smaller or same as original
      } else { 
        # if 'i' in C[i,j] is a bigger age band -> weighted average of the contact rates of participant groups
        group_span=model_agegroups$wpp_agegroup_low[i_row]:model_agegroups$wpp_agegroup_high[i_row]
        agegroup_weights=orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span])
        # if 'j' is same/smaller -> contact rate with original group proportionally divided
        if (model_agegroups$wpp_agegroup_low[j_col]==model_agegroups$wpp_agegroup_high[j_col]) {
          f_dur=model_agegroups$duration[j_col]/orig_age_groups_duration[model_agegroups$wpp_agegroup_high[j_col]]
          C_m[i_row,j_col]=
            sum((orig_age_groups_sizes[group_span]/sum(orig_age_groups_sizes[group_span]))*C_m_full[group_span,j_col])*f_dur
          # print("'i' group (C[i,j]) is bigger age band than original & 
          #       'j' is same as original group"); print(c(i_row,j_col))
          
        } else {
          # if 'j' larger -> SUM of contacts by contact groups
          j_group_span=model_agegroups$wpp_agegroup_low[j_col]:model_agegroups$wpp_agegroup_high[j_col]
          # contact numbers BY participant groups we need to AVERAGE (proport to pop size) over
          # contact numbers BY contact groups we need to SUM
          C_m[i_row,j_col]=sum(agegroup_weights*unlist(lapply(group_span, function(x) sum(C_m_full[x,j_group_span]))))
          # print("'i' group (C[i,j]) is bigger age band than original & 
          #       'j' is bigger than original group"); print(c(i_row,j_col))
          
        }
      }

      } 
    }
  C_m 
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# funcn create reciprocal matrix ----------------------------------------------------------
fun_recipr_contmatr <- function(C_m_full,age_group_sizes){
  all_perms=merge(1:nrow(C_m_full),1:nrow(C_m_full)) # permutations(n=nrow(C_m_full),r=2,repeats.allowed=T)
  N_tot=sum(age_group_sizes)
  C_m_full_symm=matrix(0,nrow=nrow(C_m_full),ncol=nrow(C_m_full))
  for (k in 1:nrow(all_perms)) { 
    i=all_perms[k,1]; j=all_perms[k,2]
    C_m_full_symm[i,j]=(C_m_full[i,j] + C_m_full[j,i]*(age_group_sizes[j]/age_group_sizes[i]))/2
  }
  colnames(C_m_full_symm)=colnames(C_m_full); rownames(C_m_full_symm)=rownames(C_m_full) 
  C_m_full_symm
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# age structure of country ----------------------------------------------------------
fun_cntr_agestr <- function(i_cntr,i_year,age_low_vals,age_high_vals){
  age_groups=data.frame(age_low=seq(0,75,5), age_high=c(seq(4,74,5),100))
  if (!any((.packages()) %in% "wpp2019")) {library(wpp2019)}; if (!exists("popF")) {data("pop")}
  cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],
                         values=popF[popF$name %in% i_cntr,i_year] + popM[popM$name %in% i_cntr,i_year])
  agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
  N_tot=cntr_agestr$values[agegr_truthvals]
  N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals])
  N_tot=N_tot*1e3; # N_tot
  data.frame(age_low=age_low_vals, age_high=age_high_vals,values=N_tot, 
             duration=(age_high_vals-age_low_vals)+1) %>% mutate(proportion=values/sum(values))
}

