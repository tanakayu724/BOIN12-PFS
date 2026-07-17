####################################################################
# BOIN12-PFS operating characteristic simulation
# This program implements the BOIN12-PFS design.
####################################################################

setwd("")
source("boin12_pfs_functions.R")



BOIN12_PFS_sim <- function(
    ndose, toxm, effm, pfsm,
    ncohort1, cohortsize,
    pi, beta0, betae, betat, lambda, kappa,
    gammad,
    etfl1, addfl1
) {
  
  # BOIN12-PFS dose-finding simulation
  #
  # ndose      : number of dose levels
  # toxm       : target toxicity probability
  # effm       : target efficacy probability
  # pfsm       : target 6-month PFS rate
  # ncohort1   : maximum number of cohorts
  # cohortsize : number of patients per cohort
  # pi         : joint toxicity-efficacy probability matrix
  # beta0      : baseline log-hazard for progression
  # betae      : efficacy effect on progression hazard
  # betat      : toxicity effect on progression hazard
  # lambda     : Weibull baseline hazard parameter
  # kappa      : Weibull shape parameter
  # gammad     : dose-specific log-hazard effect

  # etfl1      : DLT accessment window (months)
  # addfl1     : additional follow-up after enrollment completion (months)
  #
  # Returns:
  #   data.stage  : patient-level trial dataset
  #   tn          : number of patients treated at each dose
  #   t1          : number of toxicity events at each dose
  #   mtd         : estimated maximum tolerated dose
  #   data.ut     : mean utility by dose
  #   setdose     : admissible dose set
  #   stopst1     : stopping status
  #   currenttime : trial duration
  #   obd         : selected optimal biological dose
  
  
  #--------------------------------------------------
  # Utility weights and prior parameters
  #--------------------------------------------------
  bpff <- 70.6
  btox <- 17.6
  beff <- 11.8
  
  a0in <- 2
  b0in <- 0.606
  
  # pre-specified sample size cutoffs for dose finding rules
  nc1 <- 9
  nc2 <- 12
  nc3 <- 6
  nc4 <- 9
  
  #--------------------------------------------------
  # BOIN toxicity boundaries
  #--------------------------------------------------
  tx0 <- toxm
  
  tx1 <- toxm * 1.4
  boundt1 <- log((1 - tx1) / (1 - tx0)) /
    log((tx0 * (1 - tx1)) / (tx1 * (1 - tx0)))
  
  tx1 <- toxm * 0.6
  boundt2 <- log((1 - tx1) / (1 - tx0)) /
    log((tx0 * (1 - tx1)) / (tx1 * (1 - tx0)))
  
  #--------------------------------------------------
  # Initialize trial status
  #--------------------------------------------------
  cdose <- 1
  setdose <- 1:ndose
  
  tn <- rep(0, ndose)
  t1 <- rep(0, ndose)
  
  data.stage1 <- data.frame()
  
  stopst1 <- 0
  currenttime <- 0
  
  entimelist <- make_entime(months = 100, average_per_month = 3) # Not used. However, removing it changes the random-number stream.
  
  #--------------------------------------------------
  # Dose finding simulation
  #--------------------------------------------------
  for (i in 1:ncohort1) {
    
    # Generate patient-level outcomes for the current cohort
    data.temp <- Gendata.c_wib(
      cohortsize = cohortsize,
      d          = cdose,
      pi         = pi,
      lambda     = lambda,
      kappa      = kappa,
      beta0      = beta0,
      betae      = betae,
      betat      = betat,
      gammad     = gammad
    )
    
    data.temp$cohortnum <- i
    
    # Update toxicity counts
    tn[cdose] <- tn[cdose] + cohortsize
    t1[cdose] <- t1[cdose] + sum(data.temp$tox)
    
    pt <- t1 / tn
    
    # Generate enrollment times
    data.temp$entime <- currenttime +
      make_entime(months = 10, average_per_month = 3)[1:cohortsize]
    
    data.temp$aentime <- ifelse(
      data.temp$entime < currenttime,
      currenttime,
      data.temp$entime
    )
    
    currenttime <- max(data.temp$aentime) + etfl1
    
    data.stage1 <- bind_rows(data.stage1, data.temp)
    data.stage1$obspreiod <- currenttime - data.stage1$aentime
    
    # Calculate interim ORR and utility quantities
    listpe <- cal_tite_orr(
      data  = data.stage1,
      ndose = ndose
    )
    
    listpfs <- calut_pfs_orr(
      data  = data.stage1,
      ndose = ndose,
      bpff  = bpff,
      btox  = btox,
      beff  = beff,
      toxm  = toxm,
      effm  = effm,
      pfsm  = pfsm,
      a0in  = a0in,
      b0in  = b0in
    )
    
    previous_dose <- cdose
    
    #------------------------------------------------
    # Safety and efficacy admissibility
    #------------------------------------------------
    str_t <- rep(1, ndose)
    str_e <- rep(1, ndose)
    
    for (d in 1:ndose) {
      str_t[d] <- pbeta(
        toxm,
        1 + t1[d],
        1 + tn[d] - t1[d]
      )
      
      str_e[d] <- 1 - pbeta(
        effm,
        1 + listpe$e1[d],
        1 + tn[d] - listpe$e1[d]
      )
    }
    
    for (j in 1:ndose) {
      if (str_t[j] < 0.1) {
        setdose <- setdose[setdose < j]
      }
      
      if (str_e[j] < 0.1) {
        setdose <- setdose[setdose != j]
      }
    }
    
    if (length(setdose) == 0) {
      stopst1 <- 1
      break
    }
    
    #------------------------------------------------
    # Dose-finding rule
    #------------------------------------------------
    if (pt[cdose] >= boundt1) {
      
      cdose <- changedose(
        cdose   = cdose,
        setdose = setdose,
        dein    = "down"
      )
      
    } else {
      
      if (
        cdose < max(setdose) &&
        tn[cdose + 1] <= 0 &&
        tn[cdose] >= nc1
      ) {
        
        cdose <- changedose(
          cdose   = cdose,
          setdose = setdose,
          dein    = "up"
        )
        
      } else if (
        pt[cdose] <= boundt2 &&
        cdose < max(setdose) &&
        tn[changedose(cdose = cdose, setdose = setdose, dein = "up")] <= nc3 &&
        tn[cdose] >= nc2
      ) {
        
        cdose <- changedose(
          cdose   = cdose,
          setdose = setdose,
          dein    = "up"
        )
        
      } else {
        
        if (pt[cdose] <= boundt2 || tn[cdose] <= nc4) {
          
          candidate_doses <- c(
            changedose(cdose = cdose, setdose = setdose, dein = "down"),
            cdose,
            changedose(cdose = cdose, setdose = setdose, dein = "up")
          )
          
        } else {
          
          candidate_doses <- c(
            changedose(cdose = cdose, setdose = setdose, dein = "down"),
            cdose
          )
        }
        
        candidate_doses <- sort(unique(candidate_doses), na.last = TRUE)
        
        best_doses <- candidate_doses[
          listpfs$prorut[candidate_doses] ==
            max(listpfs$prorut[candidate_doses])
        ]
        
        if (length(best_doses) > 1) {
          cdose <- sample(best_doses, 1)
        } else {
          cdose <- best_doses
        }
      }
    }
    
    # Stop if the same dose has accumulated enough patients
    if (tn[previous_dose] >= 21 && previous_dose == cdose) {
      break
    }
  }
  
  #--------------------------------------------------
  # Final follow-up update
  #--------------------------------------------------
  currenttime <- currenttime + addfl1
  data.stage1$obspreiod <- currenttime - data.stage1$aentime
  
  listpp6 <- calc_pfs_landmark(
    data.stage1 = data.stage1,
    pptime      = 6,
    ndose       = ndose
  )
  
  listpe <- cal_tite_orr(
    data  = data.stage1,
    ndose = ndose
  )
  
  listut <- calut_pfs_orr(
    data  = data.stage1,
    ndose = ndose,
    bpff  = bpff,
    btox  = btox,
    beff  = beff,
    toxm  = toxm,
    effm  = effm,
    pfsm  = pfsm,
    a0in  = a0in,
    b0in  = b0in
  )
  
  ut <- listut$ut
  
  #--------------------------------------------------
  # Final admissible set
  #--------------------------------------------------
  if (length(setdose) > 0) {
    
    str_t <- rep(1, ndose)
    str_e <- rep(1, ndose)
    str_p <- rep(1, ndose)
    
    for (d in 1:ndose) {
      str_t[d] <- pbeta(
        toxm,
        1 + t1[d],
        1 + tn[d] - t1[d]
      )
      
      str_e[d] <- 1 - pbeta(
        effm,
        1 + listpe$e1[d],
        1 + tn[d] - listpe$e1[d]
      )
      
      str_p[d] <- 1 - pbeta(
        pfsm,
        1 + listpp6$ess[d] - listpp6$pd[d],
        1 + listpp6$pd[d]
      )
    }
    
    for (j in 1:ndose) {
      if (str_t[j] < 0.1) {
        setdose <- setdose[setdose < j]
      }
      
      if (str_e[j] < 0.1) {
        setdose <- setdose[setdose != j]
      }
      
      if (str_p[j] < 0.1) {
        setdose <- setdose[setdose != j]
      }
    }
    
    iso <- Iso::pava(t1[tn != 0] / tn[tn != 0])
    iso <- iso + seq_along(iso) * 1e-10
    
    mtd <- min(which(abs(iso - toxm) == min(abs(iso - toxm))))
    
    setdose <- setdose[setdose <= mtd]
    setdose <- intersect(setdose, which(tn >= 9))
  }
  
  #--------------------------------------------------
  # Select optimal biological dose
  #--------------------------------------------------
  if (length(setdose) > 0) {
    
    # which(ut == max(ut[setdose]))[1],
    obd <- setdose[which.max(ut[setdose])]
    
  } else {
    
    obd <- 99
    mtd <- 99
    
    if (stopst1 == 0) {
      stopst1 <- 3
    }
  }
  
  return(
    list(
      data.stage  = data.stage1,
      tn          = tn,
      t1          = t1,
      mtd         = mtd,
      data.ut     = ut,
      setdose     = setdose,
      stopst1     = stopst1,
      currenttime = currenttime,
      obd         = obd
    )
  )
}



start_time=Sys.time()
BOIN12_PFS_OC_sim <- function(
    ttox,
    teff,
    gammad
){
  
  # Run operating characteristic simulations for the BOIN12-PFS design.
  #
  # Repeatedly simulates complete BOIN12-PFS trials and summarizes:
  #   - OBD selection percentage
  #   - Mean patient allocation by dose
  #   - Mean total sample size
  #   - Mean trial duration
  #
  # ttox       : true toxicity probabilities by dose
  # teff       : true efficacy probability matrix
  # gammad     : dose-specific log-hazard effects
  #
  # Fixed design settings:
  #   ndose      : number of dose levels
  #   cohortsize : number of patients per cohort
  #   ncohort1   : maximum number of cohorts
  #   toxm       : target toxicity probability
  #   effm       : target efficacy probability
  #   pfsm       : target 6-month PFS rate
  #   ntrial     : number of simulated trials
  #
  # Returns:
  #   dose1-doseN  : OBD selection percentage or mean patient allocation
  #   No Selection : percentage of simulation trials with no selected dose
  #   Metric       : summary metric ("% Selected" or "N Patients")
  #   Total N      : mean total sample size
  #   Duration     : mean trial duration (months)

  
  
  #--------------------------------------------------
  # Simulation settings
  #--------------------------------------------------
  set.seed(1192)
  
  ndose      <- 5
  cohortsize <- 3
  ncohort1   <- 20
  ntrial     <- 10000
  
  toxm <- 0.35
  effm <- 0.25
  pfsm <- 0.30
  
  etfl1  <- 1
  addfl1 <- 5
  
  #--------------------------------------------------
  # Outcome generation model
  #--------------------------------------------------
  beta0 <- 1
  betae <- -0.5
  betat <- -0.1
  
  lambda <- 1
  kappa  <- 1
  
  pi <- biv(
    ttox = ttox,
    teff = teff,
    c    = 0.2
  )
  
  #--------------------------------------------------
  # Storage
  #--------------------------------------------------
  selected_dose <- integer(ntrial)
  total_n <- numeric(ntrial)
  trial_duration <- numeric(ntrial)
  
  patient_counts <- data.frame()
  
  #--------------------------------------------------
  # Simulate trials
  #--------------------------------------------------
  for(i in seq_len(ntrial)){
    
    sim <- BOIN12_PFS_sim(
      ndose      = ndose,
      toxm       = toxm,
      effm       = effm,
      pfsm       = pfsm,
      ncohort1   = ncohort1,
      cohortsize = cohortsize,
      pi         = pi,
      beta0      = beta0,
      betae      = betae,
      betat      = betat,
      lambda     = lambda,
      kappa      = kappa,
      gammad     = gammad,
      etfl1      = etfl1,
      addfl1     = addfl1
    )
    
    patient_counts <- bind_rows(
      patient_counts,
      sim$data.stage %>%
        group_by(dose) %>%
        summarise(
          record_count = n(),
          .groups = "drop"
        )
    )
    
    selected_dose[i] <- sim$obd
    total_n[i] <- sum(sim$tn)
    trial_duration[i] <- sim$currenttime
  }
  
  #--------------------------------------------------
  # OBD selection percentage
  #--------------------------------------------------
  selection_summary <-
    data.frame(dose = selected_dose) %>%
    group_by(dose) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      percentage = round(
        100 * count / sum(count),
        2
      )
    ) %>%
    select(-count)
  
  selection_summary$dose <-
    paste0("dose", selection_summary$dose)
  
  selection_table <-
    selection_summary %>%
    pivot_wider(
      names_from  = dose,
      values_from = percentage
    )
  
  selection_table$Metric <- "% Selected"
  
  #--------------------------------------------------
  # Average patient allocation
  #--------------------------------------------------
  patient_summary <-
    patient_counts %>%
    group_by(dose) %>%
    summarise(
      mean_n = round(
        sum(record_count) / ntrial,
        1
      ),
      .groups = "drop"
    )
  
  patient_summary$dose <-
    paste0("dose", patient_summary$dose)
  
  enrollment_table <-
    patient_summary %>%
    pivot_wider(
      names_from  = dose,
      values_from = mean_n
    )
  
  enrollment_table$Metric <- "N Patients"
  enrollment_table$`Total N` <- round(mean(total_n), 1)
  enrollment_table$Duration <- round(mean(trial_duration), 1)
  
  #--------------------------------------------------
  # Combine output
  #--------------------------------------------------
  output <- bind_rows(
    selection_table,
    enrollment_table
  )
  
  required_cols <- c(
    paste0("dose", 1:ndose),
    "dose99",
    "Metric",
    "Total N",
    "Duration"
  )
  
  for(col in setdiff(required_cols, names(output))){
    output[[col]] <- NA
  }
  
  output <- output[, required_cols]
  
  names(output)[names(output) == "dose99"] <-
    "No Selection"
  
  output <- data.frame(
    lapply(output, as.character),
    stringsAsFactors = FALSE
  )
  
  print(output)
  
  return(output)
}

#--------------------------------------------------
# Example scenarios
#--------------------------------------------------

scenario1 <- BOIN12_PFS_OC_sim(
  ttox = c(0.07, 0.15, 0.20, 0.25, 0.45),
  teff = matrix(c(
      0.20,    0.30,    0.40,    0.50,    0.60,
      0.50,    0.45625, 0.39375, 0.33125, 0.26875,
      0.30,    0.24375, 0.20625, 0.16875, 0.13125
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.007,-2.441, -2.695,-2.963,-3.254))

scenario2 <- BOIN12_PFS_OC_sim(
  ttox = c(0.10, 0.15, 0.18, 0.20, 0.22),
  teff = matrix(c(
      0.15,    0.30,    0.55,    0.55,    0.55,
      0.55,    0.45625, 0.24375, 0.28125, 0.31875,
      0.30,    0.24375, 0.20625, 0.16875, 0.13125
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.040,-2.441,-2.599,-2.937,-3.308))


scenario3 <- BOIN12_PFS_OC_sim(
  ttox = c(0.02, 0.04, 0.06, 0.08, 0.10),
  teff = matrix(c(
      0.20,    0.30,    0.40,    0.50,    0.60,
      0.50,    0.49375, 0.45,    0.35,    0.25,
      0.30,    0.20625, 0.15,    0.15,    0.15
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.013,-2.774,-3.221,-3.160,-3.096))

scenario4 <- BOIN12_PFS_OC_sim(
  ttox = c(0.05, 0.15, 0.30, 0.40, 0.50),
  teff = matrix(c(
      0.20,    0.40,    0.40,    0.40,    0.40,
      0.50,    0.43125, 0.43125, 0.43125, 0.43125,
      0.30,    0.16875, 0.16875, 0.16875, 0.16875
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.009,-3.034,-3.019,-3.009,-2.999))



scenario5 <- BOIN12_PFS_OC_sim(
  ttox = c(0.05, 0.10, 0.20, 0.40, 0.50),
  teff = matrix(c(
      0.30,    0.40,    0.30,    0.20,    0.10,
      0.475,   0.46875, 0.49375, 0.55625, 0.60,
      0.225,   0.13125, 0.20625, 0.24375, 0.30
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.613,-3.407,-2.757,-2.481,-2.035))


scenario6 <- BOIN12_PFS_OC_sim(
  ttox = c(0.10, 0.15, 0.22, 0.30, 0.50),
  teff = matrix(c(
      0.20,    0.40,    0.60,    0.50,    0.40,
      0.51875, 0.35625, 0.19375, 0.33125, 0.46875,
      0.28125, 0.24375, 0.20625, 0.16875, 0.13125
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.183,-2.375,-2.562,-2.958,-3.367))

scenario7 <- BOIN12_PFS_OC_sim(
  ttox = c(0.07, 0.15, 0.25, 0.35, 0.50),
  teff = matrix(c(
      0.20,    0.30,    0.40,    0.50,    0.60,
      0.51875, 0.47500, 0.45000, 0.31250, 0.13750,
      0.28125, 0.22500, 0.15000, 0.18750, 0.26250
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.186,-2.602,-3.202,-2.781,-2.033))


scenario8 <- BOIN12_PFS_OC_sim(
  ttox = c(0.05, 0.08, 0.10, 0.15, 0.20),
  teff = matrix(c(
      0.30,    0.40,    0.50,    0.60,    0.30,
      0.40,    0.39375, 0.35,    0.19375, 0.45625,
      0.30,    0.20625, 0.15,    0.20625, 0.24375
    ),nrow = 3,byrow = TRUE),
  gammad = c(-1.935,-2.707,-3.158,-2.569,-2.436))

scenario9 <- BOIN12_PFS_OC_sim(
  ttox = c(0.05, 0.10, 0.15, 0.20, 0.30),
  teff = matrix(c(
      0.4000, 0.5000, 0.4000, 0.3500, 0.2000,
      0.3375, 0.2750, 0.4125, 0.5000, 0.6125,
      0.2625, 0.2250, 0.1875, 0.1500, 0.1875
    ),nrow = 3,byrow = TRUE),
  gammad = c(-2.219,-2.477,-2.865,-3.236,-2.968))

scenario10 <- BOIN12_PFS_OC_sim(
  ttox = c(0.25, 0.35, 0.55, 0.60, 0.65),
  teff = matrix(c(
      0.1000, 0.3500, 0.4000, 0.5000, 0.5000,
      0.5625, 0.3500, 0.31875, 0.2375, 0.2375,
      0.3375, 0.3000, 0.28125, 0.2625, 0.2625
    ),nrow = 3,byrow = TRUE),
  gammad = c(-1.633,-1.865,-1.993,-2.092,-2.088))

scenario11 <- BOIN12_PFS_OC_sim(
  ttox = c(0.25, 0.35, 0.55, 0.60, 0.65),
  teff = matrix(c(
      0.1000, 0.1500, 0.20000, 0.2000, 0.2500,
      0.5625, 0.5500, 0.51875, 0.5375, 0.4875,
      0.3375, 0.3000, 0.28125, 0.2625, 0.2625
    ),nrow = 3,byrow = TRUE),
  gammad = c(-1.633,-2.013,-2.136,-2.299,-2.260))


