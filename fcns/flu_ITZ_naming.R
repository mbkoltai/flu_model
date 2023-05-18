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

