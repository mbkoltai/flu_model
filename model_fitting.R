### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# FITTING

# individual simulations in fluEvidenceSynthesis
ag <- stratify_by_age(demography, limits=c(65)) # c( 43670500, 8262600 )
population <- stratify_by_risk(ag, matrix(c(0.01,0.4),nrow=1), labels = c("LowRisk", "HighRisk")) 
# c( 43233795, 4957560, 436705, 3305040)
initial.infected <- stratify_by_risk( age_groups = c(1000,1000), matrix(c(0.01,0.4),nrow=1)) # c(990, 600, 10, 400)
# contact matrix
# Polymod data is subdivided in seven age groups
poly <- polymod_uk[,c(1,2,3,9)]
poly[,3] <- rowSums(polymod_uk[,3:8])
contacts <- fluEvidenceSynthesis::contact_matrix(polymod_data=as.matrix(poly), demography=demography,
                                                 age_group_limits = c(65))
susceptibility <- c( 0.7, 0.3 ) # Different for different ages
transmissibility <- 0.17 # Same for all ages
infection_delays <- c( 0.8, 1.8 ) # 0.8 and 1.8 day.

# SIMULATE
odes <- infectionODEs(population, initial.infected, vaccine_calendar, contact_matrix=contacts, 
                      susceptibility, transmissibility, infection_delays, 7)
head(odes %>% mutate_if(is.numeric, round) )
# % infected
fraction.infected <- odes %>%
  gather(Group, Incidence, -Time) %>%
  mutate(fraction = Incidence/population[Group])

# have ContMatrix as dataframe
cm_polymod_format <- as.matrix( data.frame(list_contact_matr[[sel_cntr]]) %>%  
                                  mutate(Age=rownames(list_contact_matr[[sel_cntr]]),Weekend=0) %>%
                                  pivot_longer(!c(Age,Weekend)) %>% 
                                  mutate(name=paste0("[",gsub(",99","+",gsub("\\.",",",gsub("X","",name))),")"),
                                         Age=case_when(Age %in% "0-5" ~ 2.5,
                                                       Age %in% "5-20" ~ 12.5,
                                                       Age %in% "20-65" ~ 42.5,
                                                       Age %in% "65-99" ~ 75)) %>%
                                  pivot_wider(names_from = name,values_from = value) )

yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))

# contact matrix data from FluEvSynth
data("polymod_uk")
polymod_FluEv <- polymod_uk[,c(1,2)]
polymod_FluEv[,3] <- rowSums(polymod_uk[,c(3,4,5,6,7,8)])
polymod_FluEv[,4] <- polymod_uk[,9]

initial_parameters <- c(0.1, 0.1, 1e-5, 0.16, 0.5, 0.5, -0.15)
names(initial_parameters) <- c("epsilon_1", "epsilon_2", "psi", "transmissibility",
                               "susceptibility_1", "susceptibility_2", "initial_infected")

# is this similar to our aggreg CM from [Prem 2021]?
polymod_uk_aggreg <- polymod_uk %>% select(!Weekend) %>%
  pivot_longer(!c(Age),names_to="contact_Age_group") %>% # ,Weekend
  mutate(Age_group=factor(case_when(
    Age <5 ~ age_group_names[1],
    (Age >=5&Age<20) ~ age_group_names[2],
    (Age >=20&Age<65) ~ age_group_names[3],
    Age >=65 ~ age_group_names[4]),levels=age_group_names)) %>% 
  group_by(Age_group,contact_Age_group) %>%
  summarise(value=mean(value)) %>% # group_by(contact_Age_group) %>% 
  mutate(contact_ag_broad=
           factor(age_group_names[findInterval(as.numeric(factor(contact_Age_group,
                                                                 levels=colnames(polymod_uk)[3:9])),vec=c(3,5,7))+1],levels=age_group_names) ) %>%
  group_by(Age_group,contact_ag_broad) %>% summarise(value=sum(value)) %>% 
  arrange(Age_group,contact_ag_broad) %>% pivot_wider(names_from=contact_ag_broad,values_from=value) %>%
  column_to_rownames(var="Age_group") 
# vs 'our' contact matrix
round(list_contact_matr$`United Kingdom`,2)

inference.results <- inference(demography = demography,
                               vaccine_calendar = vaccine_calendar,
                               polymod_data = as.matrix(polymod_FluEv),
                               ili = ili_df$ili,
                               mon_pop = ili_df$total.monitored,
                               n_pos = confirmed.samples_df$positive,
                               n_samples = confirmed.samples_df$total.samples,
                               initial = initial_parameters,
                               age_groups = c(65),
                               nbatch = 1000,
                               nburn = 1000, blen = 5)
