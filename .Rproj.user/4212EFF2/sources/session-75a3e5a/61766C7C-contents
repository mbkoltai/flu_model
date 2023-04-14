# flu project
# FluNet data

# load libraries, functions
load("fcns.R")

flunet_data <- read_csv("data/VIW_FNT.csv") 
# flunet_data$`COUNTRY/AREA/TERRITORY`[1:11]
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
# time courses, data source types by diff colors
my_colors <- c("NOTDEFINED"="#F8766D","NONSENTINEL"="#00BA38","SENTINEL"="#619CFF")
pop_lim=20
for (k_reg in unique(flunet_data$WHOREGION)) {
  ff <- flunet_data %>%
    filter(COUNTRY_AREA_TERRITORY %in% (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY & 
             ISO_YEAR>=2008 & ISO_YEAR<2020 & WHOREGION %in% k_reg) %>% # ORIGIN_SOURCE
    group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
    group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK,ORIGIN_SOURCE) %>% mutate(mean_all=mean(INF_ALL,na.rm=T)) %>% ungroup() %>%
    group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% 
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(max_cntr=max(INF_ALL,na.rm=T)) 
  # plot  
  ff %>% ggplot() +
    geom_line(aes(x=ISO_WEEK,y=INF_ALL,color=factor(ORIGIN_SOURCE),group=interaction(ISO_YEAR,ORIGIN_SOURCE)),alpha=1/2) + 
    geom_line(aes(x=ISO_WEEK,y=mean_all,color=factor(ORIGIN_SOURCE),group=ORIGIN_SOURCE), size=1.2) +
    facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + labs(color="data source") +
    geom_text(data=ff %>% select(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE,max_cntr,n_year) %>% unique(),
              show.legend = FALSE,aes(x=40+3*as.numeric(factor(ORIGIN_SOURCE)),y=0.9*max_cntr,
                        label=paste0(n_year,ifelse(as.numeric(factor(ORIGIN_SOURCE))<3,", ","")),color=ORIGIN_SOURCE)) +
    scale_color_manual(values=my_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    theme_bw() + standard_theme
  # save
  ggsave(paste0("output/incidence/source_color_code/lt_",pop_lim,"M/",k_reg,"_by_datasource",".png"),
         width=36,height=18,units="cm")
}

# time course of # positives by ISO_WEEK, years overlaid, cntrs on facets, SEPARATE DATA SOURCES
for (var_source in c("SENTINEL","NONSENTINEL","NOTDEFINED")) {
pop_lim=20; start_end_yrs=c(2008,2019)
for (k_reg in unique(flunet_data$WHOREGION)) {
ff <- flunet_data %>% 
  filter(COUNTRY_AREA_TERRITORY %in% (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY & 
           ISO_YEAR>=start_end_yrs[1] & ISO_YEAR<=start_end_yrs[2] & ISO_WEEK<=52 & 
           WHOREGION %in% k_reg & ORIGIN_SOURCE %in% var_source) %>% 
  group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR) %>% mutate(mean_all=mean(INF_ALL,na.rm=T)) %>% 
  filter(mean_all>1 & n_distinct(ISO_WEEK)>20) %>%
  group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK) %>% mutate(mean_all=mean(INF_ALL,na.rm=T)) %>%
  group_by(COUNTRY_AREA_TERRITORY) %>% mutate(n_year=n_distinct(ISO_YEAR))
# 
df_text=ff %>% ungroup %>% select(COUNTRY_AREA_TERRITORY,n_year,INF_ALL) %>% group_by(COUNTRY_AREA_TERRITORY) %>% 
  filter(INF_ALL==max(INF_ALL,na.rm = T)) %>% unique()
if (nrow(ff)>0) {
ff %>% ggplot() +
  geom_line(aes(x=ISO_WEEK,y=INF_ALL,color=ISO_YEAR,group=ISO_YEAR)) + 
  geom_line(aes(x=ISO_WEEK,y=mean_all),color="black") + # mean values
  facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + ylab("# positives (all flu)") +
  geom_text(data=df_text,aes(y=INF_ALL*0.9,label=n_year),x=50) +
  scale_color_gradientn(colours=colorRampPalette(c("orange", "red"))(12)) +
  scale_x_continuous(expand=expansion(0.01,0)) + theme_bw() + standard_theme
# save
ggsave(paste0("output/incidence/lt_",pop_lim,"M/",var_source,"/",k_reg,".png"),
       width=36,height=18,units="cm") }
}
}

# explore by strain
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
for (k_reg in unique(flunet_data$WHOREGION)) {
  ff <- flunet_data %>%
    filter(COUNTRY_AREA_TERRITORY %in% (cntr_pop_2020 %>% filter(pop_size>=pop_lim))$COUNTRY_AREA_TERRITORY & 
             ISO_YEAR>=2008 & ISO_YEAR<2020 & WHOREGION %in% k_reg) %>% 
    group_by(COUNTRY_AREA_TERRITORY,ISO_YEAR,ORIGIN_SOURCE) %>% filter(mean(INF_ALL,na.rm=T)>1) %>% ungroup() %>%
    select(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE,INF_A,INF_B) %>% 
    pivot_longer(c(INF_A,INF_B),names_to="strain") %>%
    group_by(COUNTRY_AREA_TERRITORY,ISO_WEEK,ORIGIN_SOURCE,strain) %>% 
    mutate(mean_all=mean(value,na.rm=T)) %>% ungroup() %>%
    group_by(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE) %>% mutate(n_year=n_distinct(ISO_YEAR)) %>% 
    group_by(COUNTRY_AREA_TERRITORY) %>% mutate(max_cntr=max(value,na.rm=T)) %>% 
    filter(!ORIGIN_SOURCE %in% "SENTINEL") %>% 
    arrange(COUNTRY_AREA_TERRITORY,ISO_YEAR,ISO_WEEK) %>% 
    complete(ISO_YEAR = seq(min(ISO_YEAR), max(ISO_YEAR), by = 1), 
             ISO_WEEK = seq(min(ISO_WEEK), max(ISO_WEEK), by = 1),
             ORIGIN_SOURCE=c("NOTDEFINED","NONSENTINEL"),
             strain=c("INF_A","INF_B"),fill=list(value=NA))
  # plot  
  ff %>% ggplot() +
    geom_line(aes(x=ISO_WEEK,y=value,color=factor(strain),group=interaction(strain,ISO_YEAR,ORIGIN_SOURCE)),alpha=1/2) + 
    # geom_line(aes(x=ISO_WEEK,y=mean_all,color=factor(strain),
    #               group=interaction(strain,ISO_YEAR,ORIGIN_SOURCE)),alpha=1/3,size=1.2) +
    facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + labs(color="data source") +
    # geom_text(data=ff %>% select(COUNTRY_AREA_TERRITORY,ORIGIN_SOURCE,max_cntr,n_year) %>% unique(),
    #           show.legend = FALSE,aes(x=40+3*as.numeric(factor(ORIGIN_SOURCE)),y=0.9*max_cntr,
    #         label=paste0(n_year,ifelse(as.numeric(factor(ORIGIN_SOURCE))<3,", ","")),color=ORIGIN_SOURCE)) +
    scale_color_manual(values=strain_colors) + scale_x_continuous(expand=expansion(0.01,0)) + 
    theme_bw() + standard_theme
  # save
  ggsave(paste0("output/strains/lt_",pop_lim,"M/time_courses/",k_reg,"_by_datasource",".png"),
         width=36,height=18,units="cm")
}


# # on circ coordinates
# cp <- coord_polar()# theta="y"
# cp$is_free <- function() TRUE
# 
# ff %>% ggplot() +
#   geom_line(aes(x=ISO_WEEK,y=INF_ALL,color=factor(ISO_YEAR),group=ISO_YEAR)) + 
#   geom_line(aes(x=ISO_WEEK,y=mean_all),color="black") + 
#   facet_wrap(~COUNTRY_AREA_TERRITORY,scales = "free") + cp +
#   scale_x_continuous(breaks=(0:13)*4) + theme_bw() + standard_theme

# plot indiv country
# flunet_data %>% 
#   filter(COUNTRY_AREA_TERRITORY %in% "France" & ISO_YEAR %in% c(2012:2015) ) %>%
#   #  & ISO_WEEKSTARTDATE < as.Date("2020-03-01") & ISO_WEEKSTARTDATE > as.Date("2010-01-01")
#   ggplot() + 
#   geom_line(aes(x=ISO_WEEK,y=INF_ALL,color=factor(ISO_YEAR),group=ISO_YEAR)) + 
#   # geom_line(aes(x=ISO_WEEK,y=mean_all),color="black",linewidth=1.1) + 
#   facet_wrap(~COUNTRY_AREA_TERRITORY,scales="free_y") + theme_bw() + standard_theme

# geom_boxplot(aes(x=cnt_int,color=intervention,middle=norm_median*1e2,
#                  ymin=norm_CI95_low*1e2,ymax=norm_CI95_high*1e2,
#                  lower=norm_CI50_low*1e2,upper=norm_CI50_high*1e2),
#              position=position_dodge(width=dodge_val),stat="identity",
#              show.legend=ifelse(k_plot<3,F,T)) +
#   scale_color_manual(values=c("red","blue")) +
#   facet_wrap(~plot_variable) + geom_vline(xintercept=2.5,linetype="dashed",size=0.3) + # vartype
#   theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
#   geom_text(aes(x=as.numeric(interaction(intervention,country_iso))+1/2,
#                 y=max(norm_CI95_high*1e2)*1.05,group=intervention,
#                 label=ifelse(
#         intervention!="MV",gsub("^,","",paste0(orig_burden_round,
#           ifelse(grepl("YLD non",plot_variable),"*",""))),"")),
#             position=position_dodge(width=dodge_val),size=ifelse(exists("geom_text_font_size"),
#                                                                  geom_text_font_size,5)) + 
#   labs(color="",linetype="",
#        caption=ifelse(any(grepl("costs",df_plot$plot_variable)),"*pre-intervention median value",""))
