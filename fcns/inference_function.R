### FUNCTION: custom_inference ###
custom_inference <- function(
    input_demography, 
    vaccine_calendar_list, 
    input_polymod, 
    ili = NULL, 
    mon_pop = NULL, 
    # n_pos, 
    epidemics_to_fit, # passing in all epidemics instead of n_pos 
    n_samples, 
    initial,
    mapping,
    nbatch, 
    nburn, 
    blen
) {

  age.group.limits <- c(5,20,65)  # 4 age groups used in the model
  
  prop_vacc_start <- list(
    prop_vaccine_compartments = rep(0,12),
    prop_R_vaccinated =rep(0,12), 
    prop_R = rep(0,12)
  )
  
  # Define combined log likelihood which calls llikelihood for each epidemic 
  # and generates the combined likelihood.
  
  names(initial_parameters) <- c("reporting", "transmissibility","susceptibility","initial_infected")
  
  llikelihood_all <- function(pars){
    print(paste0("Step=",global_index))
    global_index <<- global_index+1
    # if(global_index==1756){browser()}
    total_ll <- 0
    all_pars <- pars
    for (e in 1:length(epidemics_to_fit)){
      pars <- all_pars[(((e-1))*6+1):(e*6)]
      vaccine_calendar <- vaccine_calendar_list[[e]]
      n_pos <- epidemics_to_fit[[e]]$data_points
      # if("NaN" %in% pars) {browser}
      # if(global_index==11652){browser()}
      ll_epidemic <-  llikelihood (pars,n_pos,vaccine_calendar)
      # browser()
      total_ll <- total_ll + ll_epidemic
      if(total_ll=="-Inf"){break}
    }
    # browser()
    return(total_ll)
  }
  
  
  # Define the actual log likelihood function
  llikelihood <- function(pars, n_pos, vaccine_calendar) {
    # browser()
    
    contacts <- fluEvidenceSynthesis::contact_matrix(
      as.matrix(input_polymod),
      input_demography, 
      age.group.limits
    ) 
    
    age.groups <- stratify_by_age(input_demography, age.group.limits)
    
    # Population size initially infected by age
    initial.infected <- rep(10^pars[4], length(age.groups))
    pars[2] <- pars[2]/100
    
    # Run simulation
    # Note that to reduce complexity 
    # we are using the same susceptibility parameter for multiple age groups
    # browser()
    odes <- incidence_function_fit(
      demography_input = input_demography,
      parameters = pars,
      calendar_input = vaccine_calendar,
      contact_ids_sample = input_polymod, 
      contacts = contacts,
      waning_rate = 0,
      vaccination_ratio_input = prop_vacc_start,
      begin_date = vaccine_calendar$dates[1], 
      end_date = vaccine_calendar$dates[2] ,  
      year_to_run = year(vaccine_calendar$dates[1]), 
      efficacy_now =rep(0,12), 
      efficacy_next=rep(0,12),
      efficacy_next2 =rep(0,12), 
      previous_summary = NA, 
      age_groups_model = age.group.limits
      )
    # browser()
    odes <- data.table(odes)

    weekly_cases <-  odes[, sum(.SD, na.rm=TRUE), by="time" ]
    
    total_ll <- 0
    for(i in 1:length(n_pos)){
      
      if (round(as.numeric(weekly_cases[i,"V1"])) < n_pos[i]) {return(as.numeric("-Inf"))} else 
      {
        weekly_ll <-  dbinom(
          x = n_pos[i], 
          size = round(as.numeric(weekly_cases[i,"V1"])),
          prob = exp(pars[1]), 
          log = T
        )
        if(is.nan(weekly_ll)) {return(as.numeric("-Inf"))}
      }
      # browser()
      total_ll <- total_ll + weekly_ll  
    }
    # browser()
    return(total_ll)
  }
  
  ### FUNCTION: llprior ### 
  llprior <- function(pars) {
    
    # 1 is reported
    # 2 is transmission
    # 3 is susceptibility
    # 4 is initinal start
# browser()
    if (
      exp(pars[1]) < 0 || exp(pars[1]) > 1 || pars[3] < 0 || pars[3] > 1  ||  pars[4] < log(0.00001) || pars[4] > 29.5  # pars[4] < log(0.01) #pars[4] < log(0.00001) || pars[4] > 29.5  # 29.5 as 0.01% of population
    ) {return(-Inf)}
    
    # lprob <- 0
    # # prior on the R0
    # # browser()
    # R0 <- fluEvidenceSynthesis::as_R0(
    #   transmission_rate = pars[2]/100,
    #   contact_matrix = input_polymod[,3:6], # contact_matrix = contacts_matrixformat,
    #   age_groups = stratify_by_age(input_demography, limits = age.group.limits)
    # )
    # 
    # lprob <- lprob + dgamma(R0-1, shape = r0_gamma_pars[1], rate =r0_gamma_pars[2], log = T)
    # lprob <- lprob + dbeta(pars[3], shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2], log = T)
    # # prior on susceptiblity
    # return(lprob)
    return(0)
  }
  ### END OF FUNCTION: llprior ###
  
  
  # Store the contact ids used during inference
  contact.ids <- list()
  
  # Run adaptive.mcmc
  mcmc.result <- adaptive.mcmc(
    lprior = llprior, 
    llikelihood =   llikelihood_all, #llikelihood, 
    nburn = nburn, 
    initial = initial,
    nbatch = nbatch, 
    blen = blen
  )

  mcmc.result$contact.ids <- t(data.frame(contact.ids))
  mcmc.result
}
### END OF FUNCTION: custom_inference ##






### FUNCTION: infectionODEs_epidemic_yearcross2 ###
infectionODEs_epidemic_yearcross2 <- function(
    population_stratified,
    initial_infected, 
    calendar_input,
    contacts_matrixformat,
    susceptibility, 
    transmissibility, 
    infection_delays, 
    interval,
    waning_rate, 
    initial_vaccinated_prop, 
    initial_Rv_prop, 
    initial_R_prop,
    begin_date, 
    end_date,
    year_to_run,
    efficacy_now, 
    efficacy_next, 
    efficacy_next2, 
    previous_summary, ...
) {
  # browser()
  # define model timings
  t <- as.numeric(seq(begin_date, end_date, interval))
  # define age group inputs
  no_groups <- length(population_stratified)
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover 
  # risk groups
  # new_cij <- matrix(rep(0,18*18), nrow = 18)
  # for (k in 1:3) {
  #   for (l in 1:3) {
  #     lk <- (k - 1)*6 + 1
  #     ll <- (l - 1)*6 + 1
  #     new_cij[lk:(lk + 6 - 1), ll:(ll + 6 - 1)] <- contacts_matrixformat
  #   }
  # }
  new_cij <- matrix(rep(0,12*12), nrow = 12)
  for (k in 1:3) {
    for (l in 1:3) {
      lk <- (k - 1)*4 + 1
      ll <- (l - 1)*4 + 1
      new_cij[lk:(lk + 4 - 1), ll:(ll + 4 - 1)] <- contacts_matrixformat
    }
  }
  
  #Assume that all R become susceptible again at the start of each posterior
  initial_R_prop <- rep(0,no_groups)
  # specify the model
  # browser()
  mod <- gen_seeiir_ag_vacc_waning$new(
    no_groups = no_groups,
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
    calendar = matrix(calendar_input$calendar, ncol = 4*3),
    gamma1 = 2/infection_delays[1],
    gamma2 = 2/infection_delays[2], 
    num_vac_start = rep(0,no_groups) # don't need to be tracked
  )
  # browser()
  #run the model
  y <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  # calculate the cumulative values
  y <- mod$transform_variables(y)$cumI
  # Returning the differences in cumulative infections from one timestep to the other
  y <- data.frame(y[2:(nrow(y)), ] - y[1:(nrow(y) - 1), ])
  #remove the last point as it's the start of the next week
  # y <- y[1:(nrow(y)-1),] # don't think we want to remove this - this is actually the infections that happened in the prior week
  # add the dates
  y$time <- seq(begin_date,length.out = nrow(y), by = interval)
  y<-data.table(y)
  # reformat the dates
  # browser()
  return(y)
  
}
### END OF FUNCTION: infectionODEs_epidemic_yearcross2 ###





### FUNCTION: incidence_function_fit ###
incidence_function_fit <- function(
    demography_input,
    parameters,
    calendar_input,
    contact_ids_sample,
    contacts,
    waning_rate,
    vaccination_ratio_input,
    begin_date, 
    end_date,  
    year_to_run, 
    efficacy_now, 
    efficacy_next,
    efficacy_next2, 
    previous_summary, 
    age_groups_model
){
  
  risk_ratios_input <- matrix(c(rep(0,8)), ncol = 4 , byrow = T) # Not using risk groups so setting this here for now.
  
  age_group_sizes <- stratify_by_age(demography_input, age_groups_model)
  
  population_stratified <- stratify_by_risk(age_group_sizes, risk_ratios_input)
  
  initial_infected <- rep(10^parameters[4],4) #6)
  initial_infected <- stratify_by_risk(initial_infected, risk_ratios_input)
  
  susceptibility <- c(0,rep(parameters[3],3))#5))
  
  # browser()
  infectionODEs_epidemic_yearcross2(
    population_stratified = population_stratified,
    initial_infected = initial_infected,
    calendar_input = calendar_input,
    contacts_matrixformat = contacts, 
    susceptibility = susceptibility,
    transmissibility = parameters[2],
    infection_delays = infection_delays, 
    interval = 7,
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
  # browser()
}
### END OF FUNCTION: incidence_function_fit ###