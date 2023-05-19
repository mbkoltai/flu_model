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
                filter(ISO_WEEKSTARTDATE>=as.Date("2010-03-01") & ISO_WEEKSTARTDATE<=as.Date("2020-03-01")) %>%
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
if (y_log_flag) {p<-p+ scale_y_log10()}; p
# SAVE
ggsave(paste0("output/plots/epid_identif/",ifelse(n_sel>1,"ALT_CNTRS/",""),
              "low_thresh",low_thresh_val, "_up_thresh",up_thresh_val,
              "_length_lim",length_lim_val,ifelse(grepl("pos",sel_variable_val),"_byposit","_bycount"),
              ifelse(y_log_flag,"_ylog",""),".png"),
       width=36,height=24,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### 
# make a list of the epidemics, aggregated by:
# - country
# - data source (metasource)
# - strain

df_epid_threshold


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
