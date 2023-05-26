custom_inference <- function(input_demography, vaccine_calendar, input_polymod, ili = NULL, 
                             mon_pop = NULL, n_pos, n_samples, initial, mapping,
                             nbatch, nburn, blen) {
  
  current.contact.ids <- seq(1,nrow(polymod_uk))
  proposed.contact.ids <- current.contact.ids
  
  # Seven age groups used in the model
  age.group.limits <- c(2,6,12,18,60)
  
  
  # Sum all populations with a certain age into their corresponding age group
  age.group.sizes.5 <- stratify_by_age(demography, age.group.limits)
  
  # if (missing(mapping))
  #   mapping <- age_group_mapping(age.group.limits, c(5,15,45,65))
  
  prop_vacc_start <- list(prop_vaccine_compartments = rep(0,18),
                          prop_R_vaccinated =rep(0,18), 
                          prop_R = rep(0,18))
  # Define the actual log likelihood function
  llikelihood <- function( pars ) {
    
    contacts <- fluEvidenceSynthesis::contact_matrix(as.matrix(input_polymod),
                                                     input_demography, age.group.limits ) 
    
    age.groups <- stratify_by_age(input_demography, 
                                  age.group.limits )
    
    
    # Population size initially infected by age and risk group
    initial.infected <- rep( 10^pars[4], 6 ) 
    # initial.infected <- stratify_by_risk(
    #   initial.infected, risk.ratios );
    pars[2] <- pars[2]/100
    # Run simulation
    # Note that to reduce complexity 
    
    # we are using the same susceptibility parameter for multiple age groups
    
    odes <- incidence_function_fit(demography_input =popthai[,2], 
                                   parameters = pars,
                                   calendar_input= vaccine_calendar,
                                   contact_ids_sample = as.matrix(polymod.thai),
                                   waning_rate = 0,
                                   vaccination_ratio_input = prop_vacc_start,
                                   begin_date = vaccine_calendar$dates[1], 
                                   end_date = vaccine_calendar$dates[2] ,  
                                   year_to_run = year(vaccine_calendar$dates[1]), 
                                   efficacy_now =rep(0,18) , 
                                   efficacy_next=rep(0,18) ,
                                   efficacy_next2 =rep(0,18), 
                                   previous_summary = NA, 
                                   age_groups_model = age.group.limits)
    
    odes <- data.table(odes)
    # Ignore times row
    odes[,Month := month(time)]
    odes[,time := NULL]
    odes[, lapply(.SD, sum, na.rm=TRUE), by="Month" ]
    monthly_cases <-  odes[, sum(.SD, na.rm=TRUE), by="Month" ]
    total_ll <- 0
    
    for(i in 1:length(n_pos)){
      
      month_ll <-  dbinom(x = n_pos[i], size = round(as.numeric(monthly_cases[i,"V1"])),
                          prob = exp(pars[1]), log = T)
      
      
      total_ll <- total_ll + month_ll  
      
    }
    return(total_ll)
  }
  llprior <- function(pars) {
    
    # 1 is reported 
    # 2 is transmission
    # 3 is susceptibility
    # 4 is initinal start
    if (exp(pars[1]) < 0 || pars[4] < log(0.00001) || pars[4] > 29.5 ){ # 29.5 as 0.01% of population
      return(-Inf)}
    
    lprob <- 0
    # prior on the R0 
    R0 <- fluEvidenceSynthesis::as_R0(transmission_rate = pars[2]/100,
                                      contact_matrix =contacts_matrixformat, 
                                      age_groups = stratify_by_age(popthai[,2], limits = 
                                                                     age.group.limits))
    lprob <- lprob + dgamma(R0-1, shape = r0_gamma_pars[1], 
                            rate =r0_gamma_pars[2], log = T)
    lprob <- lprob + dbeta(pars[3], shape1 = sus_beta_pars[1], 
                           shape2 = sus_beta_pars[2], log = T)
    
    # prior on susceptiblity
    
    
    
    return(lprob)
  }
  
  # Store the contact ids used during inference
  contact.ids <- list()
  
  # Run adaptive.mcmc
  mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood = llikelihood, 
                               nburn = nburn, 
                               initial = initial,
                               nbatch = nbatch, blen = blen)
  
  mcmc.result$contact.ids <- t(data.frame(contact.ids))
  mcmc.result
}


infectionODEs_epidemic_yearcross2 <- function(population_stratified,
                                              initial_infected, 
                                              calendar_input,
                                              contacts_matrixformat,
                                              susceptibility, transmissibility, infection_delays, interval,waning_rate, 
                                              initial_vaccinated_prop, initial_Rv_prop, initial_R_prop,
                                              begin_date, end_date,
                                              year_to_run,
                                              efficacy_now, efficacy_next, efficacy_next2, 
                                              previous_summary, ...
) {
  
  # define model timings
  t <- as.numeric(seq(begin_date, end_date, interval))
  # define age group inputs
  no_groups <- length(population_stratified)
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover 
  # risk groups
  new_cij <- matrix(rep(0,18*18), nrow = 18)
  for (k in 1:3) {
    for (l in 1:3) {
      lk <- (k - 1)*6 + 1
      ll <- (l - 1)*6 + 1
      new_cij[lk:(lk + 6 - 1), ll:(ll + 6 - 1)] <- contacts_matrixformat
    }
  }
  
  #Assume that all R become susceptible again at the start of each posterior
  initial_R_prop <- rep(0,no_groups)
  # specify the model
  
  mod <- gen_seeiir_ag_vacc_waning$new(no_groups = no_groups,
                                       cij = new_cij,
                                       trans = transmissibility,
                                       pop = population_stratified,
                                       I0 = initial_infected,
                                       V0 = initial_vaccinated_prop,
                                       R0 = initial_R_prop,
                                       RV0 = initial_Rv_prop,
                                       susc = rep(susceptibility,3),
                                       alpha = calendar_input$efficacy[1:no_groups],
                                       omega = waning_rate,
                                       dates = calendar_input$dates,
                                       calendar = matrix(calendar_input$calendar, ncol = 6*3),
                                       gamma1 = 2/infection_delays[1],
                                       gamma2 = 2/infection_delays[2], 
                                       num_vac_start = rep(0,no_groups) # don't need to be tracked
                                       
  )
  
  #run the model
  y <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  # calculate the cumulative values
  y <- mod$transform_variables(y)$cumI
  # Returning the differences in cumulative infections from one timestep to the other
  y <- data.frame(y[2:(nrow(y)), ] - y[1:(nrow(y) - 1), ])
  #remove the last point as it's the start of the next week
  y <- y[1:(nrow(y)-1),]
  # add the dates
  y$time <- seq(begin_date,length.out = nrow(y), by =interval)
  y<-data.table(y)
  # reformat the dates
  
  return(y)
  
}

# the incidence function
incidence_function_fit <- function(demography_input, 
                                   parameters,
                                   calendar_input,
                                   contact_ids_sample,
                                   waning_rate,
                                   vaccination_ratio_input,
                                   begin_date, 
                                   end_date ,  
                                   year_to_run, 
                                   efficacy_now  , 
                                   efficacy_next ,
                                   efficacy_next2 , 
                                   previous_summary, 
                                   age_groups_model){
  
  age_group_sizes <- stratify_by_age(demography_input, age_groups_model)
  
  population_stratified <- stratify_by_risk(age_group_sizes, risk_ratios_input)
  
  initial_infected <- rep(10^parameters[4], 6)
  initial_infected <- stratify_by_risk(initial_infected, risk_ratios_input)
  
  susceptibility <- c(0,rep(parameters[3],5))
  
  infectionODEs_epidemic_yearcross2(population_stratified = population_stratified,
                                    initial_infected = initial_infected,
                                    calendar_input = calendar_input,
                                    contacts_matrixformat,
                                    susceptibility = susceptibility,
                                    transmissibility = parameters[2],
                                    infection_delays = infection_delays, interval = 7,
                                    waning_rate = waning_rate,
                                    initial_vaccinated_prop = unlist(vaccination_ratio_input[[1]]),
                                    initial_Rv_prop = unlist(vaccination_ratio_input[[2]]),
                                    initial_R_prop = unlist(vaccination_ratio_input[[3]]),
                                    begin_date = begin_date,
                                    end_date = end_date,
                                    year_to_run = year_to_run,
                                    efficacy_now = efficacy_now,
                                    efficacy_next = efficacy_next,
                                    efficacy_next2 = efficacy_next2,
                                    previous_summary = previous_summary
  )
  
}