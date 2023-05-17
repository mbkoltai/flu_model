# LOAD FLUNET DATA
flunet_data <- read_csv("data/VIW_FNT.csv") 
# POPULATION
cntr_pop_2020 = pop[-(1:29),] %>% select(name,`2020`) %>% rename(pop_size=`2020`) %>% 
  mutate(pop_size=round(pop_size/1e3,2)) %>% rename(COUNTRY_AREA_TERRITORY=name)

flunet_data_UK_summed <- left_join(flunet_data %>% 
                                     mutate(COUNTRY_AREA_TERRITORY=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),
                                                                          "United Kingdom",COUNTRY_AREA_TERRITORY),
         COUNTRY_CODE=ifelse(grepl("United Kingdom",COUNTRY_AREA_TERRITORY),"UK",COUNTRY_CODE)) %>%
  filter(COUNTRY_AREA_TERRITORY %in% gsub(", England","",unique(flu_ITZ_clusters$country_altern_name)) & 
                                              ISO_YEAR>=2008 & ISO_YEAR<2020) %>%
                                     group_by(COUNTRY_AREA_TERRITORY,COUNTRY_CODE,ISO_YEAR,ISO_WEEK,ORIGIN_SOURCE) %>%
                 summarise(INF_A=sum(INF_A),INF_B=sum(INF_B),SPEC_PROCESSED_NB=sum(SPEC_PROCESSED_NB)) %>% ungroup(),
                                   flu_ITZ_clusters %>% select(country,country_altern_name) %>% 
                                     rename(COUNTRY_AREA_TERRITORY=country_altern_name) %>% unique() ) %>% 
  mutate(country=ifelse(COUNTRY_AREA_TERRITORY %in% "United Kingdom","United Kingdom",country))

####


flu_ITZ_clusters <- read_csv("data/flu_ITZ_clusters.csv") %>% 
  pivot_longer(cols=c(kmeans_cluster_name,hi_cluster_ward_name),
               names_to="method",values_to="cluster_name") %>%
  mutate(cluster_number=ifelse(grepl("kmeans",method),kmeans_cluster_num,hi_cluster_ward_num),
         method=gsub("_name","",method)) %>% 
  select(!c(kmeans_cluster_num,hi_cluster_ward_num))

# naming inconsistencies
incons_names <- unique(flu_ITZ_clusters$country)[!unique(flu_ITZ_clusters$country) %in% 
                                                   unique(flunet_data$COUNTRY_AREA_TERRITORY)]
part_match <- unlist(lapply(incons_names, 
                            function(x) which(grepl(x,unique(flunet_data$COUNTRY_AREA_TERRITORY)))))
# unique(flunet_data$COUNTRY_AREA_TERRITORY)[part_match]

flu_ITZ_clusters <- flu_ITZ_clusters %>% 
  mutate(country_altern_name=country,
         country_altern_name=case_when(
           country %in% "Bolivia" ~ "Bolivia (Plurinational State of)",
           country %in% "Venezuela" ~ "Venezuela (Bolivarian Republic of)",
           country %in% "United States" ~ "United States of America",
           country %in% "Russia" ~ "Russian Federation",
           country %in% "United Kingdom" ~ "United Kingdom, England",
           country %in% "Iran" ~ "Iran (Islamic Republic of)",
           country %in% "Cote d'Ivoire" ~ "Côte d'Ivoire",
           country %in% "Turkey" ~ "Türkiye",
           country %in% "Laos" ~ "Lao People's Democratic Republic",
           .default=country)) # TRUE ~ country
