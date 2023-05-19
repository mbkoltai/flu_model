# source functions, call libraries
source("fcns/fcns.R")
# load data
source("fcns/load_flunet.R")

# FluEvidenceSynthesis
library(fluEvidenceSynthesis)

# load filtered data
n_sel=1
df_posit_sel_cntrs <- read_csv(c("output/df_positivity_counts_sel_cntrs.csv",
                                 "output/df_positivity_counts_sel_cntrs_ALTERN.csv")[n_sel])

# IDENTIFY THE EPIDEMICS
length_lim_val=8; low_thresh_val=0.5; up_thresh_val=0.75
sel_variable_val=c("value","positivity")[2]
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

# there are no epidemics identified based on sentinel data only, 
# so lets use the season limits from the nonsentinel (this includes nonsentinel and notdefined in FluNet)
# because this has higher counts, more years etc

# df_epid_lims %>% filter(metasource %in% "NONSENT")

# lets first do fitting for 1 country, 1 strain, 1 datatype
data_fitting <- df_epid_threshold %>% 
  filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT") %>%
  select(!c(flu_peak,over_peak,flu_included,over_inclusion,seq))
if (all(data_fitting$epidem_inclusion==(data_fitting$epid_index>0))) {
  data_fitting <- data_fitting %>% select(!epidem_inclusion)}

# df_epid_lims %>% filter(country %in% "Canada" & STRAIN %in% "INF_A" & metasource %in% "NONSENT")

### ### ### ### ### ### ### ### ### ### ### ### 
# obtain contact matrices and aggregate to our age groups:
# [0-5), [5-18), [18-65), [65-]

# load CONTACT MATRICES from [Prem 2021]
# https://github.com/kieshaprem/synthetic-contact-matrices/tree/master/output/syntheticmatrices
load("data/contact_all.rdata"); # setdiff(ls(), existing_objects)
# matrices as `contact_all$GHA`
df_cntr_table = data.frame(country=unique(df_posit_sel_cntrs$country),
                          COUNTRY_CODE=unique(df_posit_sel_cntrs$COUNTRY_CODE)) %>%
                mutate(country_sub=case_when(
                  COUNTRY_CODE %in% "UK" ~ "GBR",
                  COUNTRY_CODE %in% "AUS" ~ "NZL", TRUE ~ COUNTRY_CODE))
# extract matrices from `contact_all`
cm_list <- lapply(df_cntr_table$country_sub, function(x) contact_all[[x]])
names(cm_list) <- df_cntr_table$country

# standard_age_groups <- fun_cntr_agestr(country_sel,i_year="2020",seq(0,75,5),c(seq(4,74,5),99))
# popul_struct=fcn_cntr_fullpop(n_year="2020",country_sel)
# # RSV age groups (population data from wpp2019)
# rsv_age_groups <- fun_rsv_agegroups(standard_age_groups,popul_struct,
#                                     rsv_age_groups_low=c(0,0.5,1,1.5, 2,3,4, 5,15, 45, 65),
#                                     rsv_age_group_sizes=c(rep(0.4,4),rep(0.9,3), 9, 29, 19, 34))
# # rsv_age_groups$value=rsv_age_groups$value*67e6/sum(rsv_age_groups$value)
# ons_2020_midyear_estimates_uk <- read_csv(here::here("repo_data/ons_2020_midyear_estimates_uk.csv")) %>% 
#   mutate(age_num=as.numeric(gsub("\\+","",age)))
# low_inds<-findInterval(rsv_age_groups$age_low,ons_2020_midyear_estimates_uk$age_num)
# high_inds <- findInterval(rsv_age_groups$age_low+rsv_age_groups$duration-0.1,
#                           ons_2020_midyear_estimates_uk$age_num)
# rsv_age_groups$value <- unlist(lapply(1:length(low_inds), 
#         function(x) sum(ons_2020_midyear_estimates_uk$value[low_inds[x]:high_inds[x]])*ifelse(rsv_age_groups$duration[x]<1,
#                                                         rsv_age_groups$duration[x],1) ))

# modify contact matrix to correspond to our age groups
C_m_merged_nonrecipr=fun_create_red_C_m(C_m_full=C_m_polymod,
              model_agegroups=model_age_groups,
              orig_age_groups_duration=standard_age_groups$duration,
              orig_age_groups_sizes=standard_age_groups$values)
# make it reciprocal for the larger group
C_m=fun_recipr_contmatr(C_m_full = C_m_merged_nonrecipr,
                        age_group_sizes = model_age_groups$stationary_popul)

# only the UK is in polymod
# list_CMs <- lapply(array(sel_cntrs_per_cluster), function(x)
#                    contact_matrix(polymod, countries=x, age.limits = c(0, 5, 18, 65))$matrix,
#                    symmetric=TRUE)
# there are CMs for China and Thailand
data.frame(list_surveys()) %>% filter(grepl(paste0(array(sel_cntrs_per_cluster),collapse="|"),title))

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
