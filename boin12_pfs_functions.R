####################################################################
# Functions for BOIN12-PFS
# This script contains supporting functions for the BOIN12-PFS
#
# Required packages:
#   mvtnorm, dplyr, ReIns, tidyverse
####################################################################

library(mvtnorm)  
library(dplyr)    
library(ReIns)    
library(tidyverse)

biv=function(ttox,teff,c){
  # Generate joint probabilities of toxicity and efficacy outcomes
  # using a Gaussian copula approach.
  #
  # ttox : vector of DLT probabilities at each dose level
  # teff : 3 × ndose matrix of efficacy probabilities
  #        Row 1 = PR/CR probability
  #        Row 2 = SD probability
  #        Row 3 = PD probability
  # c    : correlation coefficient between toxicity and efficacy
  #
  # Returns:
  # A 6 × ndose matrix containing joint probabilities
  #
  # Row 1 : P(Tox = 0, Eff = PD)
  # Row 2 : P(Tox = 0, Eff = SD)
  # Row 3 : P(Tox = 0, Eff = PR/CR)
  # Row 4 : P(Tox = 1, Eff = PD)
  # Row 5 : P(Tox = 1, Eff = SD)
  # Row 6 : P(Tox = 1, Eff = PR/CR)
  #
  # Columns correspond to dose levels.
  
  at <- qnorm(ttox)
  ae <- qnorm(teff[1,])
  lambda <- ae-qnorm(teff[1,]+teff[2,])
  re <- matrix(0,nrow=6,ncol=length(ttox))
  cor <- matrix( c(1,c,c,1),nrow=2     )
  for(i in 1: length(ttox)){
    re[1,i]<-pmvnorm(  upper=c(0,lambda[i]), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=0,eff=PD
    re[2,i]<-pmvnorm( lower=c( -Inf, lambda[i]  ),  upper=c(0,0), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=0,eff=SD
    re[3,i]<-pmvnorm( lower=c( -Inf, 0  ),  upper=c(0,Inf), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=0,eff=CR/PR
    re[4,i]<-pmvnorm( lower=c(0, -Inf),  upper=c(Inf,lambda[i]), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=1, eff=PD
    re[5,i]<-pmvnorm( lower=c( 0, lambda[i]  ),  upper=c(Inf,0), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=1, eff=SD
    re[6,i]<-pmvnorm( lower=c( 0, 0  ), mean=c( at[i], ae[i]   ), sigma=cor       )[1] ## tox=1,eff=CR/PR
  }
  return(re)
}

changedose=function(cdose,setdose,dein){
  # Dose transition function
  #
  # cdose   : current dose level
  # setdose : admissible dose vector
  # dein    : dose transition decision
  #           "up"   = escalate to next higher dose
  #           "down" = de-escalate to next lower dose
  #           "stay" = remain at current dose (or move to the nearest lower admissible dose
  #                    if cdose is not admissible)
  #
  # return  : next dose level
  
  next_dose<-99999
  if(dein=="up"){
    if(cdose<max(setdose)){
      next_dose<-min(setdose[setdose>cdose])}
    else if(cdose>=max(setdose)){
      next_dose<-max(setdose)}
  }   
  else if(dein=="down"){
    if(cdose>min(setdose)){
      next_dose<-max(setdose[setdose<cdose])} 
    else if(cdose<min(setdose)){next_dose<-cdose}
    else if(cdose==min(setdose)){{next_dose<-cdose}}
  }
  else if(dein=="stay"){
    if(cdose %in% setdose){next_dose<-cdose}
    else if(!(cdose %in% setdose)){
      if(cdose>min(setdose)){next_dose<-max(setdose[setdose<cdose])}
      else if(cdose<=min(setdose)){next_dose<-min(setdose[setdose>=cdose])}
    } 
  }
  return(next_dose)  
}



make_entime=function(months,average_per_month){
  # Generate patient enrollment times under a Poisson accrual model.
  #
  # months            : total accrual period (months)
  # average_per_month : average number of enrolled patients per month
  #
  # Returns:
  # A sorted vector of enrollment times (in months).
  # The number of enrolled patients in each month follows a
  # Poisson distribution with mean = average_per_month, and
  # enrollment times are assumed to be uniformly distributed
  # within each month.
  
  data <- data.frame()
  for (month in 1:months) {
    num_people <- rpois(1, lambda = average_per_month)
    registration_times <- month - 1 + runif(num_people, min = 0, max = 1)
    if (num_people > 0) {
      new_entries <- data.frame(
        Month = month, 
        entime=registration_times
      )
      data <- rbind(data, new_entries)
    }
  }
  data=data[order(data$entime), ]$entime
  return(data)
}



calc_pfs_landmark <- function(data.stage1, pptime, ndose){
  # Calculate landmark PFS statistics using the TITE concept
  #
  # data.stage1 : patient-level trial dataset
  # pptime      : landmark PFS time point (months)
  # ndose       : number of dose levels
  #
  # Returns:
  #   data.stage1 : updated dataset including:
  #                 pd1      = progression indicator at landmark
  #                 pending  = pending landmark outcome indicator
  #                 pd_nona  = progression indicator with missing
  #                            values replaced by 0
  #                 ess      = effective sample size contribution
  #   ess         : total effective sample size at each dose
  #   pd          : total number of progression events at each dose

  
  data.stage1$pd1 <- with(
    data.stage1,
    ifelse(obspreiod >= pptime,
      ifelse(pfs <= pptime, 1, 0),
      ifelse(pfs <= obspreiod, 1, NA)))
  data.stage1$pending <- ifelse(is.na(data.stage1$pd1), 1, 0)
  data.stage1$pd_nona <- ifelse(is.na(data.stage1$pd1), 0, data.stage1$pd1)
  data.stage1$ess <- ifelse(
    is.na(data.stage1$pd1),
    data.stage1$obspreiod / pptime,
    1
  )
  ess <- rep(0, ndose)
  pd  <- rep(0, ndose)
  ess_sum <- aggregate(ess ~ dose, data = data.stage1, sum)
  pd_sum  <- aggregate(pd_nona ~ dose, data = data.stage1, sum)
  ess[ess_sum$dose] <- ess_sum$ess
  pd[pd_sum$dose]   <- pd_sum$pd_nona
  return(list(
    data.stage1 = data.stage1,
    ess = ess,
    pd = pd
  ))
}




cal_tite_orr <- function(data, ndose, timeframe = 2){
  
  # Calculate TITE-adjusted ORR quantities
  #
  # data      : patient-level dataset
  # ndose     : number of dose levels
  # timeframe : efficacy assessment window (months)
  #
  # Returns:
  # e1           : observed responders by dose
  # e1_eventess  : TITE-adjusted responder contribution for each patient
  
  data$pd1 <- with(
    data,
    ifelse(
      obspreiod >= 12,
      ifelse(pfs <= 12, 1, 0),
      ifelse(pfs <= obspreiod, 1, NA)
    )
  )
  
  data$e1 <- with(
    data,
    ifelse(
      obspreiod >= timeframe,
      ifelse(efftime <= timeframe, 1, 0),
      ifelse(efftime <= obspreiod, 1, NA)
    )
  )
  
  data$e1 <- ifelse(
    !is.na(data$pd1) &
      data$pd1 == 1 &
      (is.na(data$e1) | data$e1 != 1),
    0,
    data$e1
  )

  data$e1_nona <- ifelse(is.na(data$e1), 0, data$e1)
  
  data$ess <- ifelse(
    is.na(data$e1),
    data$obspreiod / timeframe,
    1
  )

  e1  <- rep(0, ndose)
  ess <- rep(0, ndose)
  
  agg_e1  <- aggregate(e1_nona ~ dose, data = data, sum)
  agg_ess <- aggregate(ess ~ dose, data = data, sum)
  
  e1[agg_e1$dose]   <- agg_e1$e1_nona
  ess[agg_ess$dose] <- agg_ess$ess
  
  prepp <- ifelse(ess > 0, e1 / ess, 0)
  
  data$prepp <- prepp[data$dose]
  
  data$eventess <- ifelse(
    data$prepp == 1 & data$ess == 1,
    0,
    data$prepp * (1 - data$ess) /
      (1 - data$prepp * data$ess)
  )
  
  data$eventess[is.na(data$eventess)] <- 0
  
  data$e1_eventess <- data$e1_nona + data$eventess
  
  return(list(
    e1 = e1,
    e1_eventess = data$e1_eventess
  ))
}



calut_pfs_orr <- function(data, ndose,
                          bpff, btox, beff,
                          toxm, effm, pfsm,
                          a0in = 2, b0in = 0.606){
  
  
  # Calculate dose-level PFS utility and posterior probability
  # that utility exceeds the target utility.
  #
  # Args:
  #   data  : patient-level dataset
  #   ndose : number of dose levels
  #   bpff  : utility weight for PFS
  #   btox  : utility weight for absence of toxicity
  #   beff  : utility weight for efficacy (ORR)
  #   toxm  : target toxicity rate
  #   effm  : target efficacy rate
  #   pfsm  : target 6-month PFS rate
  #   a0in  : prior alpha parameter for exponential-gamma model
  #   b0in  : prior beta parameter for exponential-gamma model
  #
  # Returns:
  #   ut     : mean utility by dose
  #   prorut : posterior probability that utility exceeds target utility

  
  
  tn <- rep(0, ndose)
  tab_n <- aggregate(rep(1, nrow(data)) ~ dose, data = data, sum)
  tn[tab_n$dose] <- tab_n[, 2]
  
  listorr <- cal_tite_orr(
    data = data,
    ndose = ndose
  )
  
  listpfs <- calc_pfs_landmark(
    data.stage1 = data,
    pptime = 12,
    ndose = ndose
  )
  
  t12d <- listpfs$data.stage1
  
  t12d$ess_include <- ifelse(
    !is.na(t12d$pd1) & t12d$pd1 == 1,
    pmin(t12d$pfs / 12, 1),
    t12d$ess
  )
  
  t12d$ess_include[is.na(t12d$ess_include)] <- 0
  
  ess_include <- rep(0, ndose)
  agg_ess_include <- aggregate(
    ess_include ~ dose,
    data = t12d,
    sum
  )
  ess_include[agg_ess_include$dose] <- agg_ess_include$ess_include
  
  a0 <- a0in + listpfs$pd
  b0 <- b0in + ess_include
  
  t12d$a0 <- a0[t12d$dose]
  t12d$b0 <- b0[t12d$dose]
  
  t12d$e1_eventess <- listorr$e1_eventess
  
  t12d$e1_eventess <- ifelse(
    !is.na(t12d$pd1) & t12d$eff == 0 & t12d$pd1 == 1,
    0,
    t12d$e1_eventess
  )
  
  t12d$pfs12 <- pmin(t12d$pfs / 12, 1)
  
  t12d$utpfs <- ifelse(
    t12d$pending == 0,
    t12d$pfs12 * bpff +
      beff * t12d$e1_eventess +
      (1 - t12d$tox) * btox,
    bpff * t12d$ess +
      (1 - t12d$tox) * btox +
      beff * t12d$e1_eventess +
      ((t12d$b0 * bpff) / (t12d$a0 - 1)) *
      (1 - (t12d$b0 / (t12d$b0 + 1 - t12d$ess))^(t12d$a0 - 1))
  )
  
  t12d$utpfs <- t12d$utpfs / 100
  t12d$utpfs[is.na(t12d$utpfs)] <- 0
  
  ut_t <- rep(0, ndose)
  agg_ut <- aggregate(
    utpfs ~ dose,
    data = t12d,
    sum
  )
  ut_t[agg_ut$dose] <- agg_ut$utpfs
  
  ut <- ut_t / (tn + 0.001)
  
  pfsm_t=(1-pfsm^2)/(-2*log(pfsm))
  
  targetut <- pfsm_t * bpff +
    beff * effm +
    (1 - toxm) * btox
  
  targetutb <- targetut / 100 + (1 - targetut / 100) / 2
  
  prorut <- 1 - pbeta(
    targetutb,
    1 + ut_t,
    1 + tn - ut_t
  )
  return(list(
    ut = ut,
    prorut = prorut
  ))
}





Gendata.c_wib=function(cohortsize,d,pi,lambda,kappa,beta0,betae,betat,
                    gammad,eff_shape=2,eff_scale=2,pd_shape=2,pd_scale=2,timeframe=2){
  
  # Generate patient-level toxicity, efficacy and PFS outcomes
  #
  # cohortsize : number of patients in the cohort
  # d          : dose level
  # pi         : joint multinomial probability matrix
  # lambda     : Weibull baseline hazard parameter
  # kappa      : Weibull shape parameter
  # beta0      : baseline log-hazard
  # betae      : efficacy effect on hazard
  # betat      : toxicity effect on hazard
  # gammad     : dose-specific hazard effect
  # eff_shape  : Weibull shape parameter for efficacy time in first 2 month
  # eff_scale  : Weibull scale parameter for efficacy time in first 2 month
  # pd_shape   : Weibull shape parameter for PD time in first 2 month
  # pd_scale   : Weibull scale parameter for PD time in first 2 month
  # timeframe  : efficacy assessment window (months)
  #
  # Returns:
  #   dose      : assigned dose level
  #   tox       : toxicity indicator
  #   sd        : stable disease indicator
  #   eff       : efficacy/response indicator
  #   efftime   : response time
  #   pd        : progression indicator
  #   pfs       : progression-free survival time
  #
  

  datate.temp <- rmultinom(cohortsize, 1, pi[, d])
  
  tox.temp <- datate.temp[4, ] + datate.temp[5, ] + datate.temp[6, ]
  sd.temp  <- datate.temp[2, ] + datate.temp[5, ]
  eff.temp <- datate.temp[3, ] + datate.temp[6, ]
  pd.temp  <- datate.temp[1, ] + datate.temp[4, ]
  

  pfs.temp     <- rep(0, cohortsize)
  efftime.temp <- rep(0, cohortsize)

  
  for (l in 1:cohortsize) {
    

    if (eff.temp[l] == 1) {
      
      ha <- exp(
        beta0 +
          betae +
          betat * tox.temp[l] +
          gammad[d]
      )
      
      prog.time <- rweibull(
        1,
        shape = kappa,
        scale = (1 / (lambda * ha))^(1 / kappa)
      )
      
      pfs.temp[l] <- prog.time + timeframe

      efftime.temp[l] <- rtweibull(
        1,
        shape = eff_shape,
        scale = eff_scale,
        endpoint = timeframe
      )
      
    } else if (sd.temp[l] == 1) {
      
      ha <- exp(
        beta0 +
          betat * tox.temp[l] +
          gammad[d]
      )
      
      prog.time <- rweibull(
        1,
        shape = kappa,
        scale = (1 / (lambda * ha))^(1 / kappa)
      )
      
      pfs.temp[l] <- prog.time + timeframe

      efftime.temp[l] <- 9999
      

    } else {
      
      efftime.temp[l] <- 9999
      
      pfs.temp[l] <- rtweibull(
        1,
        shape = pd_shape,
        scale = pd_scale,
        endpoint = timeframe
      )
    }
  }
  

  dataf <- data.frame(
    dose    = d,
    tox     = tox.temp,
    sd      = sd.temp,
    eff     = eff.temp,
    efftime = efftime.temp,
    pd      = pd.temp,
    pfs     = pfs.temp
  )
  
  return(dataf)
}
