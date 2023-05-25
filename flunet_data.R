# flu project
# FluNet data

# load libraries, functions
# this file loads libraries and some custom-made functions
source("fcns/fcns.R")

# LOAD DATA
# FLUNET
flunet_data <- read_csv("data/VIW_FNT.csv") 
# POPULATION
cntr_pop_2020 = pop[-(1:29),] %>% select(name,`2020`) %>% rename(pop_size=`2020`) %>% 
  mutate(pop_size=round(pop_size/1e3,2)) %>% rename(COUNTRY_AREA_TERRITORY=name)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT #1
### This plots shows number of positive samples per year per country
# 
# number of datapoints or mean positives (by year) by COUNTRY, YEAR, ORIGIN_SOURCE
# mycolors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
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
  xlab("") + ylab(c("mean # positives","# datapoints")[k_var]) + # scale_fill_manual(values=mycolors) + 
  theme_bw() + standard_theme + geom_vline(xintercept=(2008:2019)+1/2,size=1/3) # unique(flunet_data$WHOREGION)
if (k_var==1) {p <- p + geom_point(data=dummy_df,aes(x=ISO_YEAR,y=yvar),color=NA,fill=NA)}; p
# save
ggsave(paste0("output/datasource_stats/lt_",pop_lim,"m_pop/",
              c("mean_detects","datapoint")[k_var],"/",k_reg,"_by_source_year",".png"),
         width=36,height=18,units="cm")
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


#####
# time courses

# number of specimens received vs processed
cor(flunet_data$SPEC_RECEIVED_NB, flunet_data$SPEC_PROCESSED_NB,use = "complete.obs")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT #2
# time courses, data source types by diff colors
source_colors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
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
    scale_color_manual(values=source_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    ylab(ifelse(flag_positivity,"positivity","# positives")) +  theme_bw() + standard_theme
  # save
  ggsave(paste0("output/plots/",ifelse(positivity,"positivity/","incidence/"),
                "source_color_code/lt_",pop_lim,"M/",k_reg,"_by_datasource",".png"),
         width=36,height=18,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT #3
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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Explore by STRAIN

# some calculations on correlations btwn different variables
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT #4
# Separated by strain (colour coding), also by data source (ORIGIN_SOURCE), 
# but sources are not colour-coded. Thicker lines are means per STRAIN-ORIGIN_SOURCE
strain_colors <- c("INF_A"="#F8766D","INF_B"="#00BA38")
pop_lim=20
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# these are some calculations on correlation btwn columns in FluNet data
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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# clustering: from Chen 2023 (https://doi.org/10.1016/j.ijid.2023.02.002)

flu_ITZ_clusters <- read_csv("data/flu_ITZ_clusters.csv") %>% 
  pivot_longer(cols=c(kmeans_cluster_name,hi_cluster_ward_name),
               names_to="method",values_to="cluster_name") %>%
  mutate(cluster_number=ifelse(grepl("kmeans",method),kmeans_cluster_num,hi_cluster_ward_num),
         method=gsub("_name","",method)) %>% 
  select(!c(kmeans_cluster_num,hi_cluster_ward_num))

# calling this file to correct naming inconsistencies...
source("fcns/flu_ITZ_naming.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT #5
# PLOT by STRAIN + by clusters: separated by strain and data source, each plot shows one ITZ

strain_colors <- c("INF_A (NOTDEFINED)"="#F8766D","INF_A (NONSENTINEL)"="darkred",
                   "INF_B (NOTDEFINED)"="#00BA38","INF_B (NONSENTINEL)"="darkgreen")
# plot positivity or counts (in the dataframe called `value`)
varname=c("positivity","value")[1]
# these are the 2 clustering methods used in https://doi.org/10.1016/j.ijid.2023.02.002 paper
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
# PLOT #6
# show all curves for a cluster in the same panel, all clusters on the same plot
strain_colors <- c("INF_A (NOTDEF)"="#F8766D","INF_A (NONSENT)"="darkred","INF_A (SENT)"="orange",
                   "INF_B (NOTDEF)"="#00BA38","INF_B (NONSENT)"="darkgreen","INF_B (SENT)"="turquoise")
# show positivity (not number counts)
varname=c("positivity","value")[1]

for (k_method in c("kmeans_cluster","hi_cluster_ward")){ # [1]

subset_flunet <- left_join(flunet_data %>%
  filter(COUNTRY_AREA_TERRITORY %in% flu_ITZ_clusters$country_altern_name & 
           ISO_YEAR>=2008 & ISO_YEAR<2020),
  flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
    select(country_altern_name,cluster_name) %>% 
    rename(COUNTRY_AREA_TERRITORY=country_altern_name)) %>%
  group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% 
  filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
  select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,
         INF_A,INF_B,SPEC_PROCESSED_NB,cluster_name) %>% 
  pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
  mutate(positivity=value/SPEC_PROCESSED_NB,
         STRAIN_SOURCE=paste0(STRAIN," (",gsub("DEFINED","DEF",gsub("SENTINEL","SENT",ORIGIN_SOURCE)),")")) %>%
  group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
  complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1), 
           ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),fill=list(value=NA)) %>% filter(!is.na(cluster_name))
# average over the years
mean_weekly_df <- subset_flunet %>% filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
  group_by(ISO_WEEK,STRAIN_SOURCE,cluster_name,STRAIN) %>% 
  summarise(mean_all=mean(!!sym(varname),na.rm=T),n_sample=sum(!is.na(!!sym(varname))),
            n_cntr=n_distinct(COUNTRY_AREA_TERRITORY),n_yr=n_distinct(ISO_YEAR)) %>% 
  mutate(mean_all=ifelse(n_sample>2,mean_all,NA)) %>% filter(!is.na(cluster_name))
# averages across weeks, years and countries  
mean_all_df <- mean_weekly_df %>% group_by(STRAIN_SOURCE,cluster_name) %>%
  summarise(max_cluster=max(mean_all,na.rm=T),n_sample=round(mean(n_sample)),
         n_cntr=round(mean(n_cntr)),n_yr=round(mean(n_yr)))

# plot
subset_flunet %>% filter(positivity<=1 | is.na(positivity)) %>%
  ggplot() + facet_grid(STRAIN_SOURCE~cluster_name,scales="free_y") + 
  geom_line(aes(x=ISO_WEEK,y=positivity*100,color=factor(STRAIN_SOURCE),
                group=interaction(ISO_YEAR,COUNTRY_AREA_TERRITORY,STRAIN_SOURCE)),alpha=1/6) + 
  geom_line(data=mean_weekly_df,aes(x=ISO_WEEK,y=mean_all*100,
                             group=STRAIN_SOURCE,color=factor(STRAIN_SOURCE)),linewidth=1.1) +
  geom_line(data=mean_weekly_df,aes(x=ISO_WEEK,y=mean_all*100,
                             group=STRAIN_SOURCE),linewidth=1/3,color="black") +
  geom_text(data=mean_all_df,aes(x=40,y=80,label=paste0("<n_cntr>=",n_cntr,"\n<n_yr>=",n_yr)),size=3.25) +
  labs(color="") + ylab("positivity (%)") + guides(color=guide_legend(ncol=2)) + ggtitle(k_method) +
  scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.02,0)) + 
  scale_color_manual(values=strain_colors) + 
  theme_bw() + standard_theme + theme(legend.position="top",legend.text=element_text(size=16))
# save
ggsave(paste0("output/plots/cluster/",k_method,"/",
              gsub("value","incidence",varname),"/cntrs_sep_lines.png"),
       width=36,height=24,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT #7
# show variation and min/max instead of individual curves

# only CI50?
CI_50_FLAG=T

for (k_method in c("kmeans_cluster","hi_cluster_ward")) {

  mean_var_df <- left_join(flunet_data %>% 
                             filter(COUNTRY_AREA_TERRITORY %in% flu_ITZ_clusters$country_altern_name & 
                                      ISO_YEAR>=2008 & ISO_YEAR<2020), 
                             flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
                               select(country_altern_name,cluster_name) %>% 
                               rename(COUNTRY_AREA_TERRITORY=country_altern_name)) %>%
    group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% 
    filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
    select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,
           INF_A,INF_B,SPEC_PROCESSED_NB,cluster_name) %>% 
    pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
    mutate(positivity=value/SPEC_PROCESSED_NB,
           STRAIN_SOURCE=paste0(STRAIN," (",gsub("DEFINED","DEF",gsub("SENTINEL","SENT",ORIGIN_SOURCE)),")")) %>%
    group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1),
             ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
             fill=list(value=NA)) %>% filter(!is.na(cluster_name)) %>%
  # calculate averages ...
    filter(ISO_WEEK<53 & (positivity<=1 |is.na(positivity) )) %>% 
  group_by(ISO_WEEK,STRAIN_SOURCE,cluster_name,STRAIN) %>% 
  summarise(n_cntr=n_distinct(COUNTRY_AREA_TERRITORY),n_year=n_distinct(ISO_YEAR),
            mean_all=mean(!!sym(varname),na.rm=T),median_all=median(!!sym(varname),na.rm=T),
            min_clust=min(!!sym(varname),na.rm=T),max_clust=max(!!sym(varname),na.rm=T),
            CI50_l=quantile(positivity,probs=25/100,na.rm=T),
            CI50_u=quantile(positivity,probs=75/100,na.rm=T),
            CI95_l=quantile(positivity,probs=2.5/100,na.rm=T),
            CI95_u=quantile(positivity,probs=97.5/100,na.rm=T)) %>%
  filter(!is.na(cluster_name)) %>% group_by(STRAIN_SOURCE,cluster_name) %>% 
  mutate(n_cntr_all_weeks=round(mean(n_cntr,na.rm=T)),
         max_cluster_all_weeks=max(CI95_u,na.rm=T)) %>%
  group_by(STRAIN_SOURCE) %>% mutate(max_cluster_all_weeks=max(CI50_u,na.rm=T))

# PLOT
for (mean_median_varname in c("mean_all","median_all")) {
p <- mean_var_df %>% ggplot() +
  geom_line(aes(x=ISO_WEEK,y=get(mean_median_varname)*100,
                group=STRAIN_SOURCE,color=factor(STRAIN_SOURCE)),linewidth=1.1) +
  geom_line(aes(x=ISO_WEEK,y=get(mean_median_varname)*100,
                group=STRAIN_SOURCE),linewidth=1/3,color="black") +
  geom_ribbon(aes(x=ISO_WEEK,ymin=CI50_l*100,ymax=CI50_u*100,
                  group=STRAIN_SOURCE,fill=factor(STRAIN_SOURCE)),alpha=1/3,show.legend=F) +
  facet_grid(STRAIN_SOURCE~cluster_name,scales="free_y") + labs(color="") +
  geom_text(data=mean_var_df %>% select(cluster_name,STRAIN_SOURCE,
                                        n_cntr_all_weeks,max_cluster_all_weeks) %>% unique(),
            aes(x=36,y=max_cluster_all_weeks*0.8*100,
                label=paste0("<n_cntr>=",n_cntr_all_weeks)),show.legend=FALSE) +
  scale_color_manual(values=strain_colors) + scale_fill_manual(values=strain_colors) + 
  scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.02,0)) + 
  ggtitle(k_method) + guides(color=guide_legend(ncol=2)) +
  ylab(paste0("positivity (%, ",gsub("_all","",mean_median_varname),")"))+
  theme_bw() + standard_theme + 
  theme(legend.position="top",legend.text=element_text(size=12),title=element_text(size=17))
if (!CI_50_FLAG) {
  p <- p + geom_ribbon(aes(x=ISO_WEEK,ymin=CI95_l*100,ymax=CI95_u*100,
                       group=STRAIN_SOURCE,fill=factor(STRAIN_SOURCE)),alpha=1/4,show.legend=F) }; p
# save
ggsave(paste0("output/plots/cluster/",k_method,"/",gsub("value","incidence",varname),
              "/cntrs_",gsub("_all","",mean_median_varname),"_CI",ifelse(CI_50_FLAG,"50","50_95"),".png"), 
       width=36,height=24,units="cm")
  } # mean/median
} # clustering method

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# choose 'exemplar' countries

# strain_colors <- c("INF_A (NOTDEF)"="#F8766D","INF_A (NONSENT)"="darkred","INF_A (SENT)"="orange",
#                    "INF_B (NOTDEF)"="#00BA38","INF_B (NONSENT)"="darkgreen","INF_B (SENT)"="turquoise")
source_colors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
### ### ###

flunet_data_UK_summed <- left_join(
  flunet_data %>% 
  mutate(COUNTRY_AREA_TERRITORY=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),
                                       "United Kingdom",COUNTRY_AREA_TERRITORY),
         COUNTRY_CODE=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),"UK",COUNTRY_CODE)) %>%
  filter(COUNTRY_AREA_TERRITORY %in% gsub(", England","",unique(flu_ITZ_clusters$country_altern_name)) & 
           ISO_YEAR>=2008 & ISO_WEEKSTARTDATE<as.Date("2020-04-01")) %>%
  group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE) %>%
  summarise(INF_A=sum(INF_A),INF_B=sum(INF_B),SPEC_PROCESSED_NB=sum(SPEC_PROCESSED_NB),
            ISO_WEEKSTARTDATE=unique(ISO_WEEKSTARTDATE)) %>% ungroup(),
  flu_ITZ_clusters %>% select(country,country_altern_name) %>% 
    rename(COUNTRY_AREA_TERRITORY=country_altern_name) %>% unique() ) %>% 
  mutate(country=ifelse(COUNTRY_AREA_TERRITORY %in% "United Kingdom","United Kingdom",country))
  
# PLOT #8
# calculate number of years with data per country per data source
for (k_method in c("kmeans_cluster","hi_cluster_ward")) {
# how many years of data by ITZ + average # of tests?
data_by_ITZ <-  left_join(
      flunet_data_UK_summed, 
      flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
            mutate(country_altern_name=ifelse(country_altern_name %in% "United Kingdom, England",
                                              "United Kingdom",country_altern_name)) %>%
            select(country_altern_name,cluster_name) %>% 
            rename(COUNTRY_AREA_TERRITORY=country_altern_name) ) %>%
  # JOIN
  pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>% 
  group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,STRAIN,ORIGIN_SOURCE,cluster_name) %>%
  complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1),
           ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
           fill=list(value=0,SPEC_PROCESSED_NB=0)) %>%
  # years with data
  group_by(COUNTRY_AREA_TERRITORY,STRAIN,ORIGIN_SOURCE) %>% mutate(n_yr=n_distinct(ISO_YEAR)) %>%
  # yearly averages
  group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,ISO_YEAR,ORIGIN_SOURCE,cluster_name,STRAIN) %>%
  summarise(n_wk_pos_nonzero=sum(value>0,na.rm=T),n_wk_test_nonzero=sum(SPEC_PROCESSED_NB>0,na.rm=T),
            mean_positives=mean(value,na.rm=T),mean_spec_proc=mean(SPEC_PROCESSED_NB,na.rm=T),
            n_yr=unique(n_yr)) %>%
  pivot_longer(c(n_wk_pos_nonzero,n_wk_test_nonzero,mean_positives,mean_spec_proc),
               names_to="data_metric",values_to="count") %>%
  group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,ORIGIN_SOURCE,cluster_name,STRAIN,data_metric) %>%
  summarise(n_yr=unique(n_yr),mean=mean(count),median=median(count),
            IQR_low=quantile(count,probs=25/100),IQR_high=quantile(count,probs=75/100),
            STRAIN_SOURCE=paste0(STRAIN," (",gsub("DEFINED","DEF",
                                             gsub("SENTINEL","SENT",ORIGIN_SOURCE)),")")) %>%
  unique() %>% rename(country=COUNTRY_AREA_TERRITORY) %>% mutate_if(is.numeric, round, digits=2)
data_by_ITZ <- fcn_shortercntr_names(data_by_ITZ) 

# subset for PLOT
sel_datametric <- c("n_wk_test_nonzero","mean_spec_proc","mean_positives")[2]
exempl_cntrs <- data_by_ITZ %>% filter(data_metric %in% sel_datametric & n_yr>=5 & (!STRAIN %in% "INF_B") ) %>%
  mutate(country_altern_name=gsub(" ","\n",country_altern_name)) %>%  
  group_by(STRAIN_SOURCE,cluster_name) %>% mutate(rank_median=rank(-median)) %>%
  group_by(country_altern_name) %>% mutate(rank_overall=mean(rank_median)) %>%
  group_by(cluster_name) %>% mutate(rank_overall=dense_rank(rank_overall)) %>% 
  filter(rank_overall<4 & median>20) %>% group_by(cluster_name) %>% mutate(n_c=n_distinct(country_altern_name))
# n of years w/ data
n_yr_labels <- exempl_cntrs %>% 
  group_by(cluster_name) %>% mutate(IQR=max(IQR_high)) %>%
  group_by(country_altern_name,STRAIN_SOURCE) %>% 
  summarise(n_yr=unique(n_yr),cluster_name=unique(cluster_name),IQR=unique(IQR)) %>% 
  group_by(country_altern_name) %>% 
  summarise(n_yr=paste0("[",paste0(n_yr,collapse = ","),"]"),
            IQR=unique(IQR),cluster_name=unique(cluster_name))
# plot
exempl_cntrs %>%
ggplot(aes(x=country_altern_name)) + 
  facet_wrap(~cluster_name,scales="free") +
  geom_errorbar(aes(ymin=IQR_low,ymax=IQR_high,group=ORIGIN_SOURCE),position=position_dodge(width=3/4),width=1/3) +
  geom_point(aes(y=median,color=ORIGIN_SOURCE,fill=ORIGIN_SOURCE),position=position_dodge(width=3/4),size=3) +
  xlab("") + ylab(ifelse(grepl("spec",sel_datametric),"<# specimens processed>",sel_datametric)) + 
  geom_vline(aes(xintercept=ifelse(n_c>1,(1:(n_c-1))+1/2,NA)),size=1/2) + 
  geom_text(data=n_yr_labels,aes(y=1.3*IQR,label=n_yr)) + # guides(color=guide_legend(ncol=2))
  scale_y_log10() + scale_color_manual(values=source_colors) + labs(color="",fill="") + 
  theme_bw() + standard_theme + 
  theme(legend.position="top",legend.text=element_text(size=13),strip.text=element_text(size=15),
        axis.text.x=element_text(size=14,angle=0),axis.text.y=element_text(size=13))
# save
ggsave(paste0("output/plots/cluster/summary_stats/",k_method,"_",sel_datametric,".png"),
       width=36,height=24,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT #9
# plot selected countries vs. their cluster's mean/median
n_sel <- c(1,2)[1]
sel_cntrs_per_cluster <- c("Africa"=c("Ghana","South Africa")[n_sel],
                           "Asia-Europe"=c("Turkey","Iran")[n_sel],
                           "Eastern and Southern Asia"="China", 
                           "Europe"=c("United Kingdom","Spain")[n_sel],
                           "Northern America"=c("Canada","United States")[n_sel],
                           "Oceania-Melanesia-Polynesia"="Australia","Southern America"="Argentina")

# intersect(flunet_data_UK_summed$country,sel_cntrs_per_cluster)

strain_metasource_colors <- c("INF_A (NONSENT)"="darkred","INF_A (SENTINEL)"="orange",
                              "INF_B (NONSENT)"="#00BA38","INF_B (SENTINEL)"="turquoise")

varname=c("positivity","value")[1]

for (CI_50_FLAG in c(T,F)) {
for (k_method in c("kmeans_cluster","hi_cluster_ward")) {
  
  # all data for selected countries
  df_posit_sel_cntrs <- left_join(
    flunet_data_UK_summed,
    flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
      mutate(country_altern_name=ifelse(country_altern_name %in% "United Kingdom, England",
                                        "United Kingdom",country_altern_name)) %>%
      select(country_altern_name,cluster_name) %>% 
      rename(COUNTRY_AREA_TERRITORY=country_altern_name) ) %>%
    filter(country %in% sel_cntrs_per_cluster) %>%
    mutate(metasource=ifelse(grepl("NOT|NON",ORIGIN_SOURCE),"NONSENT",ORIGIN_SOURCE)) %>%
    pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
    mutate(positivity=100*value/SPEC_PROCESSED_NB,
           STRAIN_SOURCE=paste0(STRAIN," (",metasource,")")) %>%
    group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1),
             ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
             fill=list(value=NA)) %>% 
    filter(ISO_WEEK<53 & !is.na(cluster_name) & positivity<=100)
  
  # mean/median for clusters
  cluster_mean_median <- left_join(
    flunet_data_UK_summed,
    flu_ITZ_clusters %>% filter(method %in% k_method) %>% 
                          mutate(country_altern_name=ifelse(country_altern_name %in% "United Kingdom, England",
                                           "United Kingdom",country_altern_name)) %>%
                          select(country_altern_name,cluster_name) %>% 
                          rename(COUNTRY_AREA_TERRITORY=country_altern_name) ) %>%
    mutate(metasource=ifelse(grepl("NOT|NON",ORIGIN_SOURCE),"NONSENT",ORIGIN_SOURCE)) %>%
    pivot_longer(c(INF_A,INF_B),names_to="STRAIN") %>%
    mutate(positivity=100*value/SPEC_PROCESSED_NB, STRAIN_SOURCE=paste0(STRAIN," (",metasource,")")) %>%
    group_by(COUNTRY_AREA_TERRITORY,STRAIN) %>% 
    complete(ISO_YEAR=seq(min(ISO_YEAR),max(ISO_YEAR),by=1),
             ISO_WEEK=seq(min(ISO_WEEK),max(ISO_WEEK),by=1),
             fill=list(value=NA)) %>% filter(!is.na(cluster_name)) %>%
    # calculate averages ...
    filter(ISO_WEEK<53 & (positivity<=100 |is.na(positivity) )) %>% 
    group_by(ISO_WEEK,STRAIN_SOURCE,cluster_name,STRAIN) %>% 
    summarise(n_cntr=n_distinct(COUNTRY_AREA_TERRITORY),n_year=n_distinct(ISO_YEAR),
              mean_all=mean(!!sym(varname),na.rm=T),median_all=median(!!sym(varname),na.rm=T),
              min_clust=min(!!sym(varname),na.rm=T),max_clust=max(!!sym(varname),na.rm=T),
              CI50_l=quantile(positivity,probs=25/100,na.rm=T),
              CI50_u=quantile(positivity,probs=75/100,na.rm=T),
              CI95_l=quantile(positivity,probs=2.5/100,na.rm=T),
              CI95_u=quantile(positivity,probs=97.5/100,na.rm=T)) %>%
    mutate_if(is.numeric, round, digits=2) %>%
    filter(!is.na(cluster_name)) %>% group_by(STRAIN_SOURCE,cluster_name) %>% 
    mutate(n_cntr_all_weeks=round(mean(n_cntr,na.rm=T)),
           max_cluster_all_weeks=max(CI95_u,na.rm=T)) %>%
    group_by(STRAIN_SOURCE) %>% mutate(max_cluster_all_weeks=max(CI50_u,na.rm=T))
  
  # medians of selected countries
  df_sel_cntr_median <- df_posit_sel_cntrs %>%
    group_by(COUNTRY_AREA_TERRITORY,country,STRAIN_SOURCE,ISO_WEEK,cluster_name) %>%
    summarise(median=median(positivity,na.rm=T),mean=mean(positivity,na.rm=T),
              CI50_l=quantile(positivity,probs=25/100,na.rm=T),
              CI50_u=quantile(positivity,probs=75/100,na.rm=T))
  # max per cntr
  df_cntrs_clusters <- df_posit_sel_cntrs %>% filter(STRAIN_SOURCE %in% "INF_A (NONSENT)") %>%
    group_by(STRAIN_SOURCE) %>% mutate(pos=max(positivity)) %>% 
    select(country,cluster_name,STRAIN_SOURCE,pos) %>% unique()
  # PLOT
  p <- df_posit_sel_cntrs %>% ggplot() +
      # country median
      geom_line(data=df_sel_cntr_median,aes(x=ISO_WEEK,y=median,color=factor(STRAIN_SOURCE)),linewidth=1.1) +
      # cluster median
      geom_line(data=cluster_mean_median, aes(x=ISO_WEEK,y=median_all,group=STRAIN_SOURCE),color="black") +
      facet_grid(STRAIN_SOURCE~cluster_name,scales="free_y") + 
      geom_text(data=df_cntrs_clusters,aes(x=5+nchar(country),y=ifelse(CI_50_FLAG,50,pos*0.9),
                label=country),show.legend=FALSE) +
      scale_color_manual(values=strain_metasource_colors) + scale_fill_manual(values=strain_metasource_colors) + 
      scale_x_continuous(expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.02,0)) + 
      ggtitle(k_method) + guides(color=guide_legend(ncol=2)) + 
      labs(color="") + ylab(paste0("positivity (%, median)")) +
      theme_bw() + standard_theme + 
      theme(legend.position="top",legend.text=element_text(size=12),title=element_text(size=17))
    
    if (CI_50_FLAG) {
      p <- p + 
        # country CI50
        geom_ribbon(data=df_sel_cntr_median, aes(x=ISO_WEEK,ymin=CI50_l,ymax=CI50_u,
                        group=STRAIN_SOURCE,fill=factor(STRAIN_SOURCE)),alpha=1/3,show.legend=F) +
        # cluster CI50
        geom_line(data=cluster_mean_median, aes(x=ISO_WEEK,y=CI50_l,group=STRAIN_SOURCE),color="darkgrey") +
        geom_line(data=cluster_mean_median, aes(x=ISO_WEEK,y=CI50_u,group=STRAIN_SOURCE),color="darkgrey")
    } else {
      p <- p + 
        # country indiv years
        geom_line(aes(x=ISO_WEEK,y=positivity,group=interaction(ISO_YEAR,ORIGIN_SOURCE),
                      color=factor(STRAIN_SOURCE)),alpha=1/2) }
  p
  # save
  ggsave(paste0("output/plots/cluster/summary_stats/sel_cntrs/",ifelse(n_sel>1,"altern_select/",""),k_method,
                  ifelse(CI_50_FLAG,"_CI50","_yrs_sep"),".png"), width=36,height=24,units="cm")
    # for (mean_median_varname in c("mean_all","median_all")) {  } # mean/median
}
}

# GO TO the file model_fitting.R (same folder)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# # fluID
# 
# fluID <- read_csv("data/VIW_FID_EPI.csv")
# 
# flu_ID_summ_stats <- fluID[,1:27] %>% 
#   pivot_longer(!c(WHOREGION,FLUSEASON,HEMISPHERE,ITZ,COUNTRY_CODE,COUNTRY_AREA_TERRITORY,
#   ISO_WEEKSTARTDATE,ISO_YEAR,ISO_WEEK,MMWR_WEEKSTARTDATE,MMWR_YEAR,MMWR_WEEK,ORIGIN_SOURCE,AGEGROUP_CODE)) %>%
#   group_by(WHOREGION,COUNTRY_AREA_TERRITORY,name) %>% 
#   summarise(value=sum(!is.na(value)),n_year=length(unique(ISO_YEAR))) %>% filter(value>0)
# 
# View(flu_ID_summ_stats %>% filter(grepl(paste(array(sel_cntrs_per_cluster), collapse = "|"), COUNTRY_AREA_TERRITORY)))
# 
# # view data
# fluID %>% filter(grepl(paste(c(array(sel_cntrs_per_cluster),"Türkiye"),collapse="|"), COUNTRY_AREA_TERRITORY)) %>%
#   select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ILI_CASES)
# 
# # plot ILI
# fluID %>% filter(grepl(paste(c(array(sel_cntrs_per_cluster),"Türkiye"),collapse="|"), COUNTRY_AREA_TERRITORY) & 
#                    !is.na(ILI_CASES) & !AGEGROUP_CODE %in% "UNKNOWN") %>%
#   group_by(COUNTRY_AREA_TERRITORY,AGEGROUP_CODE) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% filter(n_year>4) %>%
#   ggplot() + facet_wrap(COUNTRY_AREA_TERRITORY~AGEGROUP_CODE,scales="free_y") + 
#   geom_line(aes(x=ISO_WEEK,y=ILI_CASES,colour=ISO_YEAR,group=ISO_YEAR)) + 
#   theme_bw() + standard_theme
# 
# # plot ILI outpatient
# # fluID %>% filter(grepl(paste(c(array(sel_cntrs_per_cluster),"Türkiye"),collapse="|"), COUNTRY_AREA_TERRITORY) & 
# #                    !is.na(ILI_OUTPATIENTS) & ILI_OUTPATIENTS>0 & !AGEGROUP_CODE %in% "UNKNOWN") %>%
# #   group_by(COUNTRY_AREA_TERRITORY,AGEGROUP_CODE) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% filter(n_year>4) %>%
# #   ggplot() + facet_wrap(COUNTRY_AREA_TERRITORY~AGEGROUP_CODE,scales="free_y") + 
# #   geom_line(aes(x=ISO_WEEK,y=ILI_OUTPATIENTS,colour=ISO_YEAR,group=ISO_YEAR)) + theme_bw() + standard_theme
# 
# # plot SARI
# fluID %>% filter(grepl(paste(array(sel_cntrs_per_cluster),"Türkiye",collapse="|"), COUNTRY_AREA_TERRITORY) & 
#                    !is.na(SARI_CASES) & !AGEGROUP_CODE %in% "UNKNOWN") %>%
#   group_by(COUNTRY_AREA_TERRITORY,AGEGROUP_CODE) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% filter(n_year>4) %>%
#   ggplot(aes(x=ISO_WEEK,y=SARI_CASES,colour=ISO_YEAR,group=ISO_YEAR)) + 
#   facet_wrap(COUNTRY_AREA_TERRITORY~AGEGROUP_CODE,scales="free_y") + 
#   geom_line() + theme_bw() + standard_theme



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# dustbin of history
#
# %>% mutate(country_altern_name=country,
#        country_altern_name=case_when(
#          country %in% "Bolivia (Plurinational State of)" ~ "Bolivia",
#          country %in% "Venezuela (Bolivarian Republic of)" ~ "Venezuela",
#          country %in% "United States of America" ~ "USA",
#          country %in% "Russian Federation" ~ "Russia",
#          country %in% "United Kingdom, England" ~ "England",
#          country %in% "Iran (Islamic Republic of)" ~ "Iran",
#          country %in% "Côte d'Ivoire" ~ "Cote d'Iv",
#          country %in% "Democratic Republic of the Congo" ~ "DRC",
#          country %in% "Central African Republic" ~ "CAR",
#          country %in% "Lao People's Democratic Republic" ~ "Laos",
#          grepl("Tanzania",country) ~ "Tanzania",
#          grepl("Korea",country) ~ "S Korea",
#          TRUE ~ country))