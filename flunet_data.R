# flu project
# FluNet data

# load libraries, functions
source("fcns.R")

# LOAD DATA
# FLUNET
flunet_data <- read_csv("data/VIW_FNT.csv") 
# POPULATION
cntr_pop_2020 = pop[-(1:29),] %>% select(name,`2020`) %>% rename(pop_size=`2020`) %>% 
  mutate(pop_size=round(pop_size/1e3,2)) %>% rename(COUNTRY_AREA_TERRITORY=name)

# number of datapoints or mean positives (by year) by COUNTRY, YEAR, ORIGIN_SOURCE
# my_colors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
k_var=1; pop_lim=20
for (k_reg in unique(flunet_data$WHOREGION)) {
  plot_var=c("mean_inf_all","n_data")[k_var]
ff <- flunet_data %>%
  filter(COUNTRY_AREA_TERRITORY %in% (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY & 
           ISO_YEAR < 2020 & ISO_YEAR > 2007 & WHOREGION %in% k_reg) %>% 
  group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% 
  summarise(n_data=n(),mean_inf_all=mean(INF_ALL,na.rm=T))
dummy_df=data.frame(ISO_YEAR=2010,COUNTRY_AREA_TERRITORY=unique(ff$COUNTRY_AREA_TERRITORY),
                    yvar=mean(ff$mean_inf_all,na.rm=T),ORIGIN_SOURCE="NOTDEFINED")
p <- ff %>% ggplot(aes(x=ISO_YEAR,y=get(plot_var),fill=ORIGIN_SOURCE)) +
  geom_bar(stat="identity",position=position_dodge(preserve="single"),color="black",size=1/3,width=8/10) + 
  facet_wrap(~COUNTRY_AREA_TERRITORY,scales = c("free_y","fixed")[k_var]) + 
  scale_y_continuous(expand=expansion(mult=0.02)) + scale_x_continuous(breaks=2008:2019,expand=expansion(add=0.3)) + 
  scale_fill_manual(values=mycolors) + xlab("") + ylab(c("mean # positives","# datapoints")[k_var]) +
  theme_bw() + standard_theme + geom_vline(xintercept=(2008:2019)+1/2,size=1/3) # unique(flunet_data$WHOREGION)
if (k_var==1) {p <- p + geom_point(data=dummy_df,aes(x=ISO_YEAR,y=yvar),color=NA,fill=NA)}; p
# save
ggsave(paste0("output/datasource_stats/lt_",pop_lim,"m_pop/",
              c("mean_detects","datapoint")[k_var],"/",k_reg,"_by_source_year",".png"),
         width=36,height=18,units="cm")
}

#####
# time courses

# number of specimens received vs processed
cor(flunet_data$SPEC_RECEIVED_NB, flunet_data$SPEC_PROCESSED_NB,use = "complete.obs")


# time courses, data source types by diff colors
my_colors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
pop_lim=20; flag_positivity=T
for (k_reg in unique(flunet_data$WHOREGION)) {
  varname <- ifelse(flag_positivity,"positivity","INF_ALL")
  filt_cntrs <- (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY
  # filter flunet data
  ff <- flunet_data %>%  filter(COUNTRY_AREA_TERRITORY %in% filt_cntrs & ISO_YEAR>=2008 & ISO_YEAR<2020 &
                                WHOREGION %in% k_reg & SPEC_PROCESSED_NB>2 ) %>% 
    mutate(positivity=INF_ALL/SPEC_PROCESSED_NB) %>%
    group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE) %>% 
    mutate(n_year=n_distinct(ISO_YEAR)) %>% filter(!is.na(COUNTRY_AREA_TERRITORY)) %>% 
    group_by(COUNTRY_AREA_TERRITORY)  %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR), max(ISO_YEAR),by=1),
             ISO_WEEK=seq(min(ISO_WEEK), max(ISO_WEEK),by=1),
             ORIGIN_SOURCE=c("NOTDEFINED","NONSENTINEL","SENTINEL"),fill=list(value=NA))
    
  # average over the years
  mean_df <- ff %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
                    group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK,ORIGIN_SOURCE) %>% 
                    summarise(mean_all=mean(!!sym(varname),na.rm=T),
                              n_year=sum(!is.na(!!sym(varname)))) %>% 
                    mutate(mean_all=ifelse(n_year>2,mean_all,NA))
  # max values
  max_df <- ff %>% filter(!is.na(!!sym(varname)) & positivity<=1) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE) %>% 
    summarise(max_cntr=max(!!sym(varname),na.rm=T),n_year=unique(n_year)) %>%
    filter(!is.na(n_year)) %>% group_by(COUNTRY_AREA_TERRITORY) %>% mutate(max_cntr=max(max_cntr,na.rm=T)) %>%
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(n_source=n_distinct(ORIGIN_SOURCE))
  # plot
  ff %>% filter(positivity<=1 | is.na(positivity)) %>% ggplot() +
    geom_line(aes(x=ISO_WEEK,y=get(varname),color=factor(ORIGIN_SOURCE),
                  group=interaction(ISO_YEAR,ORIGIN_SOURCE)),alpha=1/2) + 
    geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all,color=factor(ORIGIN_SOURCE),
                               group=ORIGIN_SOURCE),linewidth=1.2) +
    facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + labs(color="data source") +
    geom_text(data=max_df,show.legend=FALSE,aes(x=40+3*as.numeric(factor(ORIGIN_SOURCE)),y=0.9*max_cntr,
                        label=paste0(n_year,ifelse(as.numeric(factor(ORIGIN_SOURCE))<3,", ","")),
                        color=factor(ORIGIN_SOURCE))) +
    scale_color_manual(values=my_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    ylab(ifelse(flag_positivity,"positivity","# positives")) +  theme_bw() + standard_theme
  # save
  ggsave(paste0("output/plots/",ifelse(positivity,"positivity/","incidence/"),
                "source_color_code/lt_",pop_lim,"M/",k_reg,"_by_datasource",".png"),
         width=36,height=18,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# time course of # positives by ISO_WEEK, years overlaid, cntrs on facets FOR ONE DATA SOURCE
flag_positivity=T; pop_lim=20; start_end_yrs=c(2008,2019)
varname <- ifelse(flag_positivity,"positivity","INF_ALL")
for (var_source in c("SENTINEL","NONSENTINEL","NOTDEFINED")) {
for (k_reg in unique(flunet_data$WHOREGION)) {
  # subset by region and data source
  ff <- flunet_data %>% 
           filter(COUNTRY_AREA_TERRITORY %in% (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY & 
           ISO_YEAR>=start_end_yrs[1] & ISO_YEAR<=start_end_yrs[2] & ISO_WEEK<=52 & 
           (WHOREGION %in% k_reg) & (ORIGIN_SOURCE %in% var_source) & SPEC_PROCESSED_NB>2)
  if (nrow(ff)>1) {
    ff <- ff %>% mutate(positivity=INF_ALL/SPEC_PROCESSED_NB) %>%
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% 
    filter(!is.na(COUNTRY_AREA_TERRITORY)) %>% group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE)  %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR), max(ISO_YEAR),by=1),
           ISO_WEEK=seq(min(ISO_WEEK), max(ISO_WEEK),by=1),fill=list(value=NA))
  # average over the years
  mean_df <- ff %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK) %>% 
    summarise(mean_all=mean(!!sym(varname),na.rm=T),n_year=sum(!is.na(!!sym(varname)))) %>% 
    mutate(mean_all=ifelse(n_year>2,mean_all,NA))
  # max values
  max_df <- ff %>% filter(!is.na(!!sym(varname)) & positivity<=1) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE) %>% 
    summarise(max_cntr=max(!!sym(varname),na.rm=T),n_year=unique(n_year)) %>%
    filter(!is.na(n_year)) %>% group_by(COUNTRY_AREA_TERRITORY) %>% 
    mutate(max_cntr=max(max_cntr,na.rm=T))

# plot
  ff %>% filter(positivity<=1 | is.na(positivity)) %>% 
ggplot() +
  geom_line(aes(x=ISO_WEEK,y=get(varname),color=ISO_YEAR,group=ISO_YEAR)) + 
  geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all),color="black",linewidth=1.15) + # mean values
  facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + ylab(varname) +
  geom_text(data=max_df,aes(y=max_cntr*0.9,label=n_year),x=50) + 
  ggtitle(paste0(k_reg," (data source: ",var_source,")")) +
  scale_color_gradientn(colours=colorRampPalette(c("orange", "red"))(12)) +
  scale_x_continuous(expand=expansion(0.01,0)) + theme_bw() + standard_theme
# save
# output/plots/incidence/INF_ALL/sources_separate/
ggsave(paste0("output/plots/",ifelse(flag_positivity,"positivity","incidence"),"/sources_separate/lt_",
              pop_lim,"M/",var_source,"/",k_reg,".png"), width=36,height=18,units="cm") }
}
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Explore by STRAIN
View(flunet_data %>% select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,contains(c("AH","SUBTYP","INF"))) %>% 
       arrange(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK))
# is INF_A + INF_B=INF_ALL? 99.99% correlated, so effectively yes
cor(flunet_data$INF_ALL,(flunet_data$INF_A+flunet_data$INF_B),use = "complete.obs")
# INF_A ~ AH... -> 98.7%
cor(flunet_data$INF_A,flunet_data$AOTHER_SUBTYPE+
    flunet_data$AH1N12009+flunet_data$AH1+flunet_data$AH3+flunet_data$AH5+flunet_data$AH7N9,use="complete.obs")
# INF_B ~ indiv B strains -> 98.9
cor(flunet_data$INF_B, flunet_data$BVIC_2DEL+flunet_data$BVIC_3DEL+flunet_data$BNOTDETERMINED+
    flunet_data$BVIC_NODEL+flunet_data$BVIC_DELUNK+flunet_data$BYAM,use="complete.obs")

# plot by strains
strain_colors <- c("INF_A"="#F8766D","INF_B"="#00BA38")
pop_lim=10
filt_countrs <- (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY
varname=c("positivity","value")[1]
for (k_reg in unique(flunet_data$WHOREGION)) {
  ff <- flunet_data %>%
    filter(COUNTRY_AREA_TERRITORY %in% filt_countrs & ISO_YEAR>=2008 & !ORIGIN_SOURCE %in% "SENTINEL" &
             ISO_YEAR<2020 & WHOREGION %in% k_reg) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
    select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,INF_A,INF_B,SPEC_PROCESSED_NB) %>% 
    pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
    mutate(positivity=value/SPEC_PROCESSED_NB) %>%
    group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
    complete(ISO_YEAR = seq(min(ISO_YEAR), max(ISO_YEAR), by=1), 
             ISO_WEEK = seq(min(ISO_WEEK), max(ISO_WEEK), by=1),
             ORIGIN_SOURCE=c("NOTDEFINED","NONSENTINEL"),fill=list(value=NA))
  # average over the years
  mean_df <- ff %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK,ORIGIN_SOURCE,STRAIN) %>% 
    summarise(mean_all=mean(!!sym(varname),na.rm=T),max_cntr=max(!!sym(varname),na.rm=T),
              n_year=sum(!is.na(!!sym(varname)))) %>% 
    mutate(mean_all=ifelse(n_year>2,mean_all,NA)) %>%
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(max_cntr=max(max_cntr,na.rm=T))
    max_df_n_years <- mean_df %>% select(COUNTRY_AREA_TERRITORY,STRAIN,ORIGIN_SOURCE,max_cntr,n_year) %>% 
      unique() %>% group_by(COUNTRY_AREA_TERRITORY,STRAIN,ORIGIN_SOURCE,max_cntr) %>%
      summarise(n_year=max(n_year))
  # plot  
  ff %>% filter(positivity<=1 | is.na(positivity)) %>%  
    ggplot() +
    geom_line(aes(x=ISO_WEEK,y=positivity,color=factor(STRAIN),
                  group=interaction(STRAIN,ISO_YEAR,ORIGIN_SOURCE)),alpha=1/2) + 
    geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all,color=factor(STRAIN),
                               group=interaction(STRAIN,ORIGIN_SOURCE)),size=1.2) +
    # alpha=1/3,
    facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + labs(color="data source") +
    geom_text(data=max_df_n_years, show.legend=FALSE,
              aes(color=factor(STRAIN), x=40+3*as.numeric(factor(ORIGIN_SOURCE)),
                  y=0.9*max_cntr+as.numeric(factor(STRAIN))*0.12,
                  label=paste0(n_year,ifelse(as.numeric(factor(ORIGIN_SOURCE))<2,", ","")))) +
    scale_color_manual(values=strain_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    ggtitle(k_reg) + theme_bw() + standard_theme
  # save
  ggsave(paste0("output/plots/",gsub("value","incidence",varname),"/strains/lt_",pop_lim,"M/",k_reg,".png"),
         width=36,height=18,units="cm")
}

# all flu samples (negatives incl.)
flunet_data %>%
  ggplot() + geom_point(aes(x=INF_ALL,y=INF_NEGATIVE,color=COUNTRY_AREA_TERRITORY),alpha=1/2) + 
  scale_x_log10() + scale_y_log10() +  
  facet_wrap(~WHOREGION,scales="free") + guides(color="none") + theme_bw() + standard_theme
#
cor(flunet_data$INF_NEGATIVE,flunet_data$INF_ALL,use="complete.obs")
# if we have # positive samples, then also # of all samples?
sum(is.na(flunet_data$INF_NEGATIVE) & is.na(flunet_data$INF_ALL))/nrow(flunet_data) # 18.6%
sum(!is.na(flunet_data$INF_NEGATIVE) & !is.na(flunet_data$INF_ALL))/nrow(flunet_data) # 15.3%
sum(is.na(flunet_data$INF_NEGATIVE) & !is.na(flunet_data$INF_ALL))/nrow(flunet_data) # 56.4%
sum(!is.na(flunet_data$INF_NEGATIVE) & is.na(flunet_data$INF_ALL))/nrow(flunet_data) # 9.7% 
# INF A
sum(!is.na(flunet_data$INF_NEGATIVE) & !is.na(flunet_data$INF_A))/nrow(flunet_data) # 56.4%
# INF B
sum(!is.na(flunet_data$INF_NEGATIVE) & !is.na(flunet_data$INF_B))/nrow(flunet_data) # 56.4%

# # on polar coordinates
# cp <- coord_polar()# theta="y"
# cp$is_free <- function() TRUE
# 
# ff %>% ggplot() +
#   geom_line(aes(x=ISO_WEEK,y=INF_ALL,color=factor(ISO_YEAR),group=ISO_YEAR)) + 
#   geom_line(aes(x=ISO_WEEK,y=mean_all),color="black") + 
#   facet_wrap(~COUNTRY_AREA_TERRITORY,scales = "free") + cp +
#   scale_x_continuous(breaks=(0:13)*4) + theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# clustering: from Chen 2023 (https://doi.org/10.1016/j.ijid.2023.02.002)

flu_ITZ_clusters <- read_csv("data/flu_ITZ_clusters.csv") %>% 
  pivot_longer(cols=c(kmeans_cluster_name,hi_cluster_ward_name),
               names_to="method",values_to="cluster_name") %>%
  mutate(cluster_number=ifelse(grepl("kmeans",method),kmeans_cluster_num,hi_cluster_ward_num),
         method=gsub("_name","",method)) %>% 
  select(!c(kmeans_cluster_num,hi_cluster_ward_num))

# naming inconsistencies
incons_names <- unique(flu_ITZ_clusters$country)[!unique(flu_ITZ_clusters$country) %in% unique(flunet_data$COUNTRY_AREA_TERRITORY)]
part_match <- unlist(lapply(incons_names, function(x) which(grepl(x,unique(flunet_data$COUNTRY_AREA_TERRITORY)))))
unique(flunet_data$COUNTRY_AREA_TERRITORY)[part_match]

flu_ITZ_clusters <- flu_ITZ_clusters %>% 
                      mutate(country_altern_name=country,
                             country_altern_name=case_when(country %in% "Bolivia" ~ "Bolivia (Plurinational State of)",
                                       country %in% "Venezuela" ~ "Venezuela (Bolivarian Republic of)",
                                       country %in% "United States" ~ "United States of America",
                                       country %in% "Russia" ~ "Russian Federation",
                                       country %in% "United Kingdom" ~ "United Kingdom, England",
                                       country %in% "Iran" ~ "Iran (Islamic Republic of)",
                                       country %in% "Cote d'Ivoire" ~ "Côte d'Ivoire",
                                       country %in% "Turkey" ~ "Türkiye",
                                       country %in% "Laos" ~ "Lao People's Democratic Republic",.default=country))

# PLOT by STRAIN + by clusters
strain_colors <- c("INF_A (NOTDEFINED)"="#F8766D","INF_A (NONSENTINEL)"="darkred",
                   "INF_B (NOTDEFINED)"="#00BA38","INF_B (NONSENTINEL)"="darkgreen")
# "INF_A (NOTDEFINED)"  "INF_B (NOTDEFINED)"  "INF_A (NONSENTINEL)" "INF_B (NONSENTINEL)"
varname=c("positivity","value")[1]
k_method=c("kmeans_cluster","hi_cluster_ward")[2]
for (k_reg in unique(flu_ITZ_clusters$cluster_name)) {
  # cntrs in cluster
  filt_countrs <- (flu_ITZ_clusters %>% filter(cluster_name %in% k_reg & method %in% k_method))$country_altern_name
  # pop_filter <- (cntr_pop_2020 %>% filter(pop_size>=5))$COUNTRY_AREA_TERRITORY
  # subset
  ff <- flunet_data %>%
    filter(COUNTRY_AREA_TERRITORY %in% filt_countrs & ISO_YEAR>=2008 & 
             !ORIGIN_SOURCE %in% "SENTINEL" & ISO_YEAR<2020) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% 
    filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
    select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,INF_A,INF_B,SPEC_PROCESSED_NB) %>% 
    pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
    mutate(positivity=value/SPEC_PROCESSED_NB,STRAIN_SOURCE=paste0(STRAIN," (",ORIGIN_SOURCE,")")) %>%
    group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1), 
             ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
             ORIGIN_SOURCE=c("NOTDEFINED","NONSENTINEL"),fill=list(value=NA))
  # average over the years
  mean_df <- ff %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
    group_by(ISO_WEEK,COUNTRY_AREA_TERRITORY,STRAIN_SOURCE) %>% 
    summarise(mean_all=mean(!!sym(varname),na.rm=T),max_cntr=max(!!sym(varname),na.rm=T),
              n_year=sum(!is.na(!!sym(varname)))) %>% 
    mutate(mean_all=ifelse(n_year>2,mean_all,NA)) %>%
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(max_cntr=max(max_cntr,na.rm=T))
  # maxima
  max_df_n_years <- mean_df %>% select(COUNTRY_AREA_TERRITORY,STRAIN_SOURCE,max_cntr,n_year) %>% 
    unique() %>% filter(!is.na(STRAIN_SOURCE)) %>% group_by(COUNTRY_AREA_TERRITORY,STRAIN_SOURCE,max_cntr) %>%
    summarise(n_year=max(n_year)) %>% 
    mutate(ORIGIN_SOURCE=ifelse(grepl("NONSENT",STRAIN_SOURCE),"NONSENTINEL","NOTDEFINED"))
  # plot
  ff %>% filter(positivity<=1 | is.na(positivity)) %>%
    ggplot() +
    geom_line(aes(x=ISO_WEEK,y=positivity,color=factor(STRAIN_SOURCE),
                  group=interaction(ISO_YEAR,STRAIN_SOURCE)),alpha=1/2) + 
    geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all,color=factor(STRAIN_SOURCE),
                               group=STRAIN_SOURCE),linewidth=1.2) + # alpha=1/3,
    facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + labs(color="data source") +
    geom_text(data=max_df_n_years, aes(x=40+3*as.numeric(factor(ORIGIN_SOURCE)),
                  y=0.9*max_cntr+ifelse(grepl("INF_A",STRAIN_SOURCE),0,1)*0.12,color=factor(STRAIN_SOURCE),
                  label=paste0(n_year,ifelse(as.numeric(factor(ORIGIN_SOURCE))<2,", ",""))), show.legend=FALSE) +
    scale_color_manual(values=strain_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    ggtitle(k_reg) + theme_bw() + standard_theme
  # save
  ggsave(paste0("output/plots/cluster/",k_method,"/",
                gsub("value","incidence",varname),"/strains_color_code/",gsub(" ","_",k_reg),".png"),
         width=36,height=18,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  
# show all curves for a cluster in the same panel, all clusters on the same plot

strain_colors <- c("INF_A (NOTDEFINED)"="#F8766D","INF_A (NONSENTINEL)"="darkred","INF_A (SENTINEL)"="orange",
                   "INF_B (NOTDEFINED)"="#00BA38","INF_B (NONSENTINEL)"="darkgreen","INF_B (SENTINEL)"="turquoise")

k_method=c("kmeans_cluster","hi_cluster_ward")[1]

subset_flunet <- left_join(flunet_data %>%
  filter(COUNTRY_AREA_TERRITORY %in% flu_ITZ_clusters$country_altern_name & ISO_YEAR>=2008 & ISO_YEAR<2020), 
  # !ORIGIN_SOURCE %in% "SENTINEL" & 
  flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
    select(country_altern_name,cluster_name) %>% rename(COUNTRY_AREA_TERRITORY=country_altern_name)) %>%
  group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% 
  filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
  select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,INF_A,INF_B,SPEC_PROCESSED_NB,cluster_name) %>% 
  pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
  mutate(positivity=value/SPEC_PROCESSED_NB,STRAIN_SOURCE=paste0(STRAIN," (",ORIGIN_SOURCE,")")) %>%
  group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
  complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1), 
           ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
           ORIGIN_SOURCE=c("NOTDEFINED","NONSENTINEL","SENTINEL"),fill=list(value=NA)) %>% filter(!is.na(cluster_name))
# average over the years
mean_df <- subset_flunet %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
  group_by(ISO_WEEK,STRAIN_SOURCE,cluster_name,STRAIN) %>% 
  summarise(mean_all=mean(!!sym(varname),na.rm=T),max_cntr=max(!!sym(varname),na.rm=T),
            n_year=sum(!is.na(!!sym(varname)))) %>% 
  mutate(mean_all=ifelse(n_year>2,mean_all,NA)) %>% filter(!is.na(cluster_name))

# plot
subset_flunet %>% filter(positivity<=1 | is.na(positivity)) %>%
  ggplot() + geom_line(aes(x=ISO_WEEK,y=positivity*100,color=factor(STRAIN_SOURCE),
                group=interaction(ISO_YEAR,COUNTRY_AREA_TERRITORY,STRAIN_SOURCE)),alpha=1/6) + 
  geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all*100,group=STRAIN_SOURCE,color=factor(STRAIN_SOURCE)),linewidth=1.1) +
  geom_line(data=mean_df,aes(x=ISO_WEEK,y=mean_all*100,group=STRAIN_SOURCE),linewidth=1/3,color="black") +
  facet_grid(STRAIN_SOURCE~cluster_name,scales="free_y") + labs(color="") +
  scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.02,0)) + ylab("positivity (%)") +
  scale_color_manual(values=strain_colors) + guides(color=guide_legend(ncol=2)) +
  ggtitle(k_method) + theme_bw() + standard_theme + theme(legend.position="top",legend.text=element_text(size=16))
# save
ggsave(paste0("output/plots/cluster/",k_method,"/",gsub("value","incidence",varname),"/cntrs_sep_lines.png"), 
       width=36,height=22,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### 
# show variation and min/max instead of individual curves

k_method=c("kmeans_cluster","hi_cluster_ward")[1]

mean_var_df <- subset_flunet %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
  group_by(ISO_WEEK,STRAIN_SOURCE,cluster_name,STRAIN) %>% 
  summarise(n_cntr=n_distinct(COUNTRY_AREA_TERRITORY),n_year=sum(!is.na(!!sym(varname)))/n_cntr,
            mean_all=mean(!!sym(varname),na.rm=T),
            min_clust=min(!!sym(varname),na.rm=T),max_clust=max(!!sym(varname),na.rm=T),
            CI50_l=quantile(positivity,probs=25/100,na.rm=T),CI50_u=quantile(positivity,probs=75/100,na.rm=T),
            CI95_l=quantile(positivity,probs=2.5/100,na.rm=T),CI95_u=quantile(positivity,probs=97.5/100,na.rm=T)) %>%
  filter(!is.na(cluster_name)) %>% group_by(STRAIN_SOURCE,cluster_name) %>% 
  mutate(n_cntr_all_weeks=round(mean(n_cntr,na.rm=T),1),max_cluster_all_weeks=max(CI95_u,na.rm=T)) %>%
  group_by(STRAIN_SOURCE) %>% mutate(max_cluster_all_weeks=max(CI50_u,na.rm=T))

# PLOT
CI_50_FLAG=T
p <- mean_var_df %>% ggplot() +
  geom_line(aes(x=ISO_WEEK,y=mean_all*100,group=STRAIN_SOURCE,color=factor(STRAIN_SOURCE)),linewidth=1.1) +
  geom_line(aes(x=ISO_WEEK,y=mean_all*100,group=STRAIN_SOURCE),linewidth=1/3,color="black") +
  geom_ribbon(aes(x=ISO_WEEK,ymin=CI50_l*100,ymax=CI50_u*100,
                  group=STRAIN_SOURCE,fill=factor(STRAIN_SOURCE)),alpha=1/3,show.legend=F) +
  facet_grid(STRAIN_SOURCE~cluster_name,scales="free_y") + labs(color="") +
  geom_text(data=mean_var_df %>% select(cluster_name,STRAIN_SOURCE,n_cntr_all_weeks,max_cluster_all_weeks) %>% 
              ungroup() %>% unique(),
            aes(x=40,y=max_cluster_all_weeks*0.8*100,label=paste0("<n>=",n_cntr_all_weeks)),show.legend=FALSE) +
  scale_color_manual(values=strain_colors) + scale_fill_manual(values=strain_colors) + ylab("positivity (%)") +
  scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.02,0)) + 
  ggtitle(k_method) + theme_bw() + standard_theme + guides(color=guide_legend(ncol=2)) +
  theme(legend.position="top",legend.text=element_text(size=12),title=element_text(size=17))
if (!CI_50_FLAG) {
  p <- p + geom_ribbon(aes(x=ISO_WEEK,ymin=CI95_l*100,ymax=CI95_u*100,
                       group=STRAIN_SOURCE,fill=factor(STRAIN_SOURCE)),alpha=1/4,show.legend=F) }
p
# save
ggsave(paste0("output/plots/cluster/",k_method,"/",gsub("value","incidence",varname),
              "/cntrs_means_CI",ifelse(CI_50_FLAG,"50","50_95"),".png"), 
       width=36,height=22,units="cm")

