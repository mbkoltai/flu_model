### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# FITTING
# we are using built-in fitting algorithm of FluEvSynthesis
# and building on the existing script of Thailand paper at 
# https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Fitting/run_inference.R

# incidence data from FluNet for selected countries are in
# the 'metasource' column is from merging the "NOTDEFINED" and "NONSENTINEL" columns,
# in most countries we only have "NOTDEFINED", and in Europe often only "NONSENTINEL",
# i am assuming both are from hospital-related lab testing
if (!exists("df_posit_sel_cntrs")) {
df_posit_sel_cntrs <- read_csv("output/df_positivity_counts_sel_cntrs.csv")
}
# incidence data WITH epidemics identified
df_epid_threshold <- read_csv(file="output/df_epid_threshold.csv")
# `df_epid_lims` contains the start and end dates of the epidemics
df_epid_lims

# 4x4 contact matrices are in `list_contact_matr` also save as RDS file in
list_contact_matr< - readRDS("output/list_contact_matr.RDS")
# this is a list with the countries as names: names(list_contact_matr)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# individual simulations in fluEvidenceSynthesis
# from https://blackedder.github.io/flu-evidence-synthesis/modelling.html
ag <- stratify_by_age(demography, limits=c(65)) # c( 43670500, 8262600 )
population <- stratify_by_risk(ag, matrix(c(0.01,0.4),nrow=1), labels = c("LowRisk", "HighRisk")) 
# c( 43233795, 4957560, 436705, 3305040)
initial.infected <- stratify_by_risk( age_groups = c(1000,1000), matrix(c(0.01,0.4),nrow=1)) # c(990, 600, 10, 400)
# contact matrix
# Polymod data is subdivided in seven age groups
fluEvSynth_poly <- polymod_uk[,c(1,2,3,9)]; fluEvSynth_poly[,3] <- rowSums(polymod_uk[,3:8])
fluEvSynth_CM <- fluEvidenceSynthesis::contact_matrix(polymod_data=as.matrix(fluEvSynth_poly), demography=demography,
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


# contact matrix data from FluEvSynth
data("polymod_uk")
polymod_FluEv <- polymod_uk[,c(1,2)]; polymod_FluEv[,3] <- rowSums(polymod_uk[,c(3,4,5,6,7,8)])
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
  mutate(contact_ag_broad=factor(age_group_names[findInterval(as.numeric(
                            factor(contact_Age_group,levels=colnames(polymod_uk)[3:9])),
                            vec=c(3,5,7))+1], levels=age_group_names) ) %>%
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# INFERENCE
# from https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Fitting/run_inference.R

# we need our contact matrices in this format: as.matrix(polymod_FluEv)

# have ContMatrix as dataframe
contacts_matrixformat <- as.matrix(
  data.frame(particip_age_group=rownames(list_contact_matr[[sel_cntr]]),
             Weekend=0,list_contact_matr[[sel_cntr]]) %>%  
        pivot_longer(!c(particip_age_group,Weekend),names_to = "contact_age") %>% rowwise() %>%
        mutate(contact_age=paste0("[",gsub(",100","+",gsub("\\.",",",gsub("X","",contact_age))),")"),
               min_part_age=as.numeric(strsplit(split="-",x=particip_age_group)[[1]][1]),
               max_part_age=as.numeric(strsplit(split="-",x=particip_age_group)[[1]][2])-1,
               Age=min_part_age) %>% group_by(particip_age_group,contact_age) %>% 
    complete(Age=seq(min_part_age, max_part_age,by=1)) %>% fill(value,Weekend) %>% ungroup() %>%
    select(!c(min_part_age,max_part_age,particip_age_group)) %>% relocate(Age,Weekend,.before =1) %>%
    pivot_wider(names_from = contact_age,values_from = value) )

yr_res_pop <- unlist(lapply(pop_age(wpp_age(sel_cntr, 2015))$population/5, function(x) rep(x,5)))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# from here i'm adjusting the code from 
# https://github.com/NaomiWaterlow/NextGen_Flu_Thai/blob/main/Fitting/run_inference.R
# but not working yet
###### things to specify #####

epidemic_to_run <- 2 # we want to fit jointly instead!
post_size <- 10000
thinning_steps <- 100
burn_in <- 100000
seed_to_use <- 55
save <- T

# lets first subset data for one country and NONSENT data only
data_fitting <- df_epid_threshold %>% 
  filter(grepl("NONSENT",metasource) & country %in% sel_cntr & STRAIN %in% "INF_A") %>%
  select(!c(over_peak,flu_included,over_inclusion,positivity,flu_peak,seq))
# we will want to do this for multiple epidemics, but for now lets select one
df_start_end <- (df_epid_lims %>% ungroup() %>% filter(index==1 & grepl("NONSENT",metasource) & STRAIN %in% "INF_A" & 
          (country %in% ifelse(grepl("King",sel_cntr),"UK",sel_cntr))))[,c("start_date","end_date")]
dates_to_run <- c(df_start_end$start_date,df_start_end$end_date)
# vacc calendar (no fitting in the vaccination)
vaccine_calendar <- as_vaccination_calendar(efficacy = rep(0,length(age_group_names)), 
                                            dates = dates_to_run,
                                            coverage = matrix(0,nrow=3, # length(dates_to_run), 
                                                              ncol = length(age_group_names)), 
                                            no_age_groups = length(age_group_names), no_risk_groups=1)
# create a list of the epidemics
epidemics_to_fit <- list(
  
  list(start = as.Date("2006-04-01"), 
       end = as.Date("2006-12-31"),
       initial_params = c(-10.8, 10, 0.6, 3, 0, 0), 
       data_points = unlist(fluAH1_epidemics[13:21,"flu"]), 
       type = "AH1N1"), )
