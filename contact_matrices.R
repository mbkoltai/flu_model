# source functions, call libraries
source("fcns/fcns.R")
# load data
source("fcns/load_flunet.R")

# load filtered data
# this csv file is saved version of `df_posit_sel_cntrs` dataframe in the flunet_data.R file
n_sel=1
df_posit_sel_cntrs <- read_csv(c("output/df_positivity_counts_sel_cntrs.csv",
                                 "output/df_positivity_counts_sel_cntrs_ALTERN.csv")[n_sel])

# PLOT #1
# IDENTIFY THE EPIDEMICS
# set thresholds (length of season, lower and upper percentiles to include as an epidemic)
length_lim_val=8; low_thresh_val=0.5; up_thresh_val=0.75
sel_variable_val=c("value","positivity")[2] # value means counts
df_epid_threshold <- fcn_identify_seasons(
  df_input = df_posit_sel_cntrs %>% 
                filter(ISO_WEEKSTARTDATE>=as.Date("2010-01-01") & ISO_WEEKSTARTDATE<=as.Date("2020-04-01")) %>%
                group_by(country,ISO_YEAR,ISO_WEEK,ISO_WEEKSTARTDATE,cluster_name,STRAIN,metasource) %>%
                summarise(SPEC_PROCESSED_NB=sum(SPEC_PROCESSED_NB),value=sum(value)) %>% ungroup(),
  sel_variable=sel_variable_val, source_varname="metasource",
  low_thresh=low_thresh_val,up_thresh = up_thresh_val,length_lim=length_lim_val)
# identify limits of epidemics
for (y_log_flag in c(T,F)) {
df_epid_lims <- fcn_find_bloc_lims(df_epid_threshold,log_flag=y_log_flag)

# plot
p <- df_epid_threshold %>%
  mutate(country=ifelse(grepl("King",country),"UK",country),
         value=ifelse(value==0,ifelse(y_log_flag,NA,value),value),
         source_epid=paste0(metasource,"_",ifelse(epidem_inclusion==1,"_ON","_OFF"))) %>%
  ggplot() + facet_grid(country~STRAIN,scales="free_y") + 
  geom_line(aes(x=ISO_WEEKSTARTDATE,y=value,group=metasource,colour=factor(source_epid))) + 
  xlab("") + ylab("# positive tests") + labs(colour="epidemic",fill="") +
  scale_color_manual(values=c("grey32","red","grey46","darkred")) + scale_fill_manual(values=c("red","darkred")) +
  scale_x_date(date_labels="%Y\n%m",breaks="6 month",expand=expansion(0.01,0))  +
  geom_rect(data=df_epid_lims, aes(xmin=start_date,xmax=end_date,fill=metasource,ymin=min_val,ymax=max_val),alpha=1/4) +
  theme_bw() + standard_theme + theme(legend.position="top",strip.text=element_text(size=13)) + 
  guides(color=guide_legend(nrow=2))
if (y_log_flag) {p<-p+ scale_y_log10(expand=expansion(0.01,0))}; p
# SAVE
ggsave(paste0("output/plots/epid_identif/",ifelse(n_sel>1,"ALT_CNTRS/",""),
              "low_thresh",low_thresh_val, "_up_thresh",up_thresh_val,
              "_length_lim",length_lim_val,ifelse(grepl("pos",sel_variable_val),"_byposit","_bycount"),
              ifelse(y_log_flag,"_ylog",""),".png"),
       width=36,height=24,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### 
# list of epidemics
#
# there are no epidemics identified based on sentinel data only, 
# so lets use the season limits from the nonsentinel (this includes nonsentinel and notdefined in FluNet)
# because this has higher counts, more years etc

# lets first subset data for 1 country, 1 strain, 1 datatype ---> fitting 
data_fitting <- df_epid_threshold %>% 
  filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT") %>%
  select(!c(flu_peak,over_peak,flu_included,over_inclusion,seq))
if (all(data_fitting$epidem_inclusion==(data_fitting$epid_index>0))) {
  data_fitting <- data_fitting %>% select(!epidem_inclusion)   }

# df_epid_lims %>% filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Get contact matrices and aggregate to our age groups:
# [0-5), [5-20), [20-65), [65-]
library(socialmixr)
# FluEvidenceSynthesis
library(fluEvidenceSynthesis)

# load CONTACT MATRICES from [Prem 2021]
# https://github.com/kieshaprem/synthetic-contact-matrices/tree/master/output/syntheticmatrices
load("data/contact_all.rdata") 
# matrices are `contact_all$GHA` ...
df_cntr_table = data.frame(country=unique(df_posit_sel_cntrs$country),
                          COUNTRY_CODE=unique(df_posit_sel_cntrs$COUNTRY_CODE)) %>%
                mutate(country_sub=case_when(
                  COUNTRY_CODE %in% "UK" ~ "GBR",
                  COUNTRY_CODE %in% "AUS" ~ "NZL", TRUE ~ COUNTRY_CODE))
# extract matrices from `contact_all`
cm_list <- lapply(df_cntr_table$country_sub, function(x) contact_all[[x]])
names(cm_list) <- df_cntr_table$country

# aggregate into our age groups
age_limits <- c(0,5,20,65); age_group_names <- paste0(age_limits,"-", c(age_limits[2:length(age_limits)],99))
# contact_matrix(age.limits = age_limits)

# with socialmixr
# orig_agegroups <- colnames(readRDS("data/UK_contact_matrix_sum.RDS"))
# xx <- as.numeric(gsub("\\+","",unlist(lapply(orig_agegroups, function(x) strsplit(x, split ="-")[[1]][1]))))
# cm_polymod_uk <- contact_matrix(polymod, countries="United Kingdom", age.limits=xx)$matrix
# cm_polymod_uk_merged <- contact_matrix(polymod, countries="United Kingdom", age.limits = c(0, 5, 20,65))$matrix

# we need popul sizes by our age groups
for (sel_cntr in df_cntr_table$country) {
# age group sizes we need (pop_age from socialmixr)
model_age_groups <- data.frame(agegroup_name=age_group_names, duration=diff(c(age_limits,120)), 
                               wpp_agegroup_low=c(1,2,5,14), wpp_agegroup_high=c(1,4,13,16),
                               popul=pop_age(wpp_age(sel_cntr, 2015), age.limit=age_limits)$population)
# age group population sizes corresponding to [Prem 2021] matrices
standard_age_groups <- fun_cntr_agestr(i_cntr = sel_cntr,i_year="2015",
                              age_low_vals = seq(0,75,5),age_high_vals = c(seq(4,74,5),120))

# modify contact matrix to correspond to our age groups
sel_cntr_code <- df_cntr_table$country_sub[df_cntr_table$country %in% sel_cntr]
C_m_merged_nonrecipr <- fun_create_red_C_m(C_m_full=list(cm_polymod_uk,contact_all[[sel_cntr_code]])[[2]],
                                        model_agegroups=model_age_groups,
                                        orig_age_groups_duration=standard_age_groups$duration,
                                        orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
if (!exists("list_contact_matr")) { list_contact_matr <- list()}
list_contact_matr[[sel_cntr]] <- fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                        age_group_sizes = model_age_groups$popul)
}

write_rds(list_contact_matr,"output/list_contact_matr.RDS")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# # Engl ILI dataset
# as.data.frame(ili$ili) %>% mutate(week=1:52) %>% pivot_longer(!week,names_to = "age") %>% 
#   ggplot() + geom_line(aes(x=week,y=value)) + facet_wrap(~age,scales="free_y") + 
#   theme_bw() + standard_theme
# 
# ### ###
# 
# cromer_samples_by_age <- qread("/home/lshmk17/Downloads/cromer_samples_by_age.qs")
# cromer_summ_stat <- cromer_samples_by_age %>% group_by(risk_group,subtype,outcome,age) %>% 
#   summarise(median=median(value),mean=mean(value),
#             ci95_l=quantile(value,probs=0.025),ci95_u=quantile(value,probs=0.975),
#             ci50_l=quantile(value,probs=0.25),ci50_u=quantile(value,probs=0.75))
# 
# cromer_summ_stat %>% filter(age<99) %>% # mutate(across(c(mean,ci95_l,ci95_u,ci50_l,ci50_u), ~ ifelse(. == 0, NA, .))) %>%
#   ggplot(aes(x=age)) + facet_grid(outcome~risk_group,scales="free_y") +
#   geom_line(aes(y=mean*100,color=subtype),size=1.2) + 
#   geom_ribbon(aes(ymin=ci50_l*100,ymax=ci50_u*100,fill=subtype),alpha=1/5) +
#   ylab("probability (%)") + # scale_y_log10() + 
#   theme_bw() + standard_theme # 
