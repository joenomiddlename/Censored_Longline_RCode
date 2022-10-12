# Simulation study demonstrating the censored hook competition method
# 3 species groups present with uncensored catch count distributions
# Target species is non-schooling, non-overdispersed.
# Same as Sim Exp 2 but the target species' abundance DECREASES
# The competitor's abundance either increases or decreases now
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)

# Flag saying whether or not the simulation is to be run, or
# pre-compiled results loaded
run_sim <- F

# Note that species 1 = 'aggressive'; 2 = 'other' group; 3 = 'target'
nspecies <- 3
# Specifies the bite time distributions of three species
# 'constant' = uniform(0,5), 'exp_decay' = trunc_exp(1), 'mixed' = 0.5*unif(0,5) + 0.5*trunc_exp(1)
bite_funs <- cbind(rep('constant',3),
                   c('exp_decay','mixed', 'constant'))
soak_time <- 5
n_hooks <- 800
# definition of stations and years are opposite to the manuscript
nstation <- 6
nyears <- 30
# higher saturation level -> more fish of all species
saturation_level <- c(1,2)
# The Target species' abundance always decreases in this experiment
# is the Competitor species' abundance increasing or decreasing?
mean_attract <- c('increase', 'decrease')
# SD of lognormal distributions which scale mean #attracted individuals
sd_log_bite <- c(0.8, 0.2, 0)
n_sim <- 10 # run 10 times
hook_sat_level <- 0.85 # true proportion at which saturation effects begin
cprop=0.85 # assumed proportion at which saturation effects begin
cprop2=0.95 # assumed proportion "" "" second model

# what does the capture probability of each species (linearly) decrease to at 100% saturation
# assume linear decrease from 85% saturation onwards
saturation_effects <- cbind(c(0, 0, 0),
                           c(0, 0.2, 0.8))
# create saturation function
sat_fun <- function(sat_effect, sat_level, hook_sat_level=0.85)
{
  val <- rep(1, length(sat_level))

  val[which(sat_level>hook_sat_level)] <-
    (1 - sat_effect*((sat_level[which(sat_level>hook_sat_level)]-hook_sat_level)/(1-hook_sat_level)))
  return(val)
}

# hook competition adjustment factor function
comp_factor_fun <- function(prop_hook, n_hook)
{
  prop <- prop_hook
  # if all hooks saturated - map to 1 hook
  prop[which(prop == 0)] <- 1 / n_hook[which(prop == 0)]
  return(-log(prop)/(1-prop))
}

# sample the unadjusted bite times for each of the 'attracted' fish for each fishing event
bite_samp <- function(bite_fun, n)
{
  if(bite_fun=='exp_decay')
  {
    val <- rexp(n)
  }
  if(bite_fun=='constant')
  {
    val <- runif(n, min = 0, max=soak_time)
  }
  if(bite_fun=='mixed')
  {
    val <- c(rexp(n),
             runif(n, min = 0, max=soak_time))[
               rbinom(n=n, size = 1, prob=0.5)+1
             ]
  }
  return(val)
}

# Create df for storing results
Results <- data.frame(
  nsim = rep(1:n_sim, each=nstation*2*2*2*2*4),
  sat_level = rep(rep(c('low','high'), times=n_sim),each=2*2*2*4*nstation),
  mean_attract = rep(rep(mean_attract, times=n_sim*2), each=2*2*4*nstation),
  sat_effects = rep(rep(c('no saturation','saturation'),times=n_sim*2*2), each=2*4*nstation),
  bite_fun = rep(rep(c('constant','mixed'), times=n_sim*2*2*2), each=4*nstation),
  model=rep(rep(c('naive','adjust','censored','censored_95'), times=n_sim*2*2*2), each=nstation),
  Bias=rep(0, times=n_sim*2*2*2*2*4*nstation),
  RMSE=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Coverage=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Rel_Bias=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Rel_RMSE=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Rel_Coverage=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Station=rep(1:nstation, times=n_sim*2*2*2*2*4),
  Prop_Sat_85=rep(0, times=n_sim*2*2*2*2*4*nstation),
  Prop_Sat_100=rep(0, times=n_sim*2*2*2*2*4*nstation))

# Create function for mapping results to correct row in df
results_mapper <- function(n,i,j,k,l,mod)
{
  return(
    which(Results$nsim == n & Results$sat_level == c('low','high')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$sat_effects == c('no saturation','saturation')[k] &
            Results$bite_fun == c('constant','mixed')[l] &
            Results$model == mod)
         )
}

if(run_sim)
{
  for(nsim in 1:n_sim)
  {
    print(paste0('iteration ',nsim,' out of ',n_sim))
    for(i in 1:length(saturation_level))
    {
      for(j in 1:length(mean_attract))
      {
        for(k in 1:dim(saturation_effects)[2])
        {
          if(mean_attract[j] == 'increase')
          {
            # Competitor increasing abundance
            # other species group constant abundance
            # Target species' abundance decreases
            mean_bite =
              cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                    rep(400,6),
                    c(220, 200, 180, 160, 140, 120)-100)
          }
          if(mean_attract[j] == 'decrease')
          {
            # Competitor decreasing abundance
            # other species group constant abundance
            # Target species' abundance decreasing
            mean_bite =
              cbind(c(280, 200, 180, 160, 140, 120)*saturation_level[i],
                    rep(400,6),
                    c(220, 200, 180, 160, 140, 120)-100)
          }
          for(l in 1:dim(bite_funs)[2])
          {
            saturation_effect <- saturation_effects[,k]
            # sample the number of each species that WOULD bite at each station for each year if hooks were available
            nbite <- data.frame(bites = rep(0, times=nspecies*nstation*nyears),
                                attracted = rep(0, times=nspecies*nstation*nyears),
                                species=rep(1:3, each=nstation*nyears),
                                station=rep(1:nstation, times=nspecies*nyears),
                                year=rep(rep(1:nyears, each=nstation), times=nspecies))
            nbite$attracted <- rpois(dim(nbite)[1],
                                     lambda = as.numeric(mean_bite[cbind(nbite$station,nbite$species)])*
                                       exp(rnorm(dim(nbite)[1], mean = 0, sd=sd_log_bite[nbite$species])))
            
            for(i2 in 1:nstation)
            {
              for(j2 in 1:nyears)
              {
                # bite times for Competitor
                bite_time_1 <- bite_samp(bite_funs[1,l],sum(nbite$attracted[nbite$species==1 &
                                                                              nbite$station==i2 &
                                                                              nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_1)>soak_time)
                {
                  bite_time_1[bite_time_1>soak_time] <-
                    bite_samp(bite_funs[1,l],sum(bite_time_1>soak_time))
                }
                # repeat for `other' species group
                bite_time_2 <- bite_samp(bite_funs[2,l],sum(nbite$attracted[nbite$species==2 &
                                                                              nbite$station==i2 &
                                                                              nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_2)>soak_time)
                {
                  bite_time_2[bite_time_2>soak_time] <-
                    bite_samp(bite_funs[2,l],sum(bite_time_2>soak_time))
                }
                
                # repeat for Target species
                bite_time_3 <- bite_samp(bite_funs[3,l],sum(nbite$attracted[nbite$species==3 &
                                                                              nbite$station==i2 &
                                                                              nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_3)>soak_time)
                {
                  bite_time_3[bite_time_3>soak_time] <-
                    bite_samp(bite_funs[3,l],sum(bite_time_3>soak_time))
                }
                
                # Now we sample the first n_hooks*n_hook_sat_level unadjusted
                if((length(bite_time_1) + length(bite_time_2) + length(bite_time_3)) <= n_hooks*hook_sat_level)
                {
                  nbite$bites[nbite$species==1 &
                                nbite$station==i2 &
                                nbite$year==j2] <- length(bite_time_1)
                  nbite$bites[nbite$species==2 &
                                nbite$station==i2 &
                                nbite$year==j2] <- length(bite_time_2)
                  nbite$bites[nbite$species==3 &
                                nbite$station==i2 &
                                nbite$year==j2] <- length(bite_time_3)
                }
                if((length(bite_time_1) + length(bite_time_2) + length(bite_time_3)) > n_hooks*hook_sat_level)
                {
                  species_ind <- c(rep(1, length(bite_time_1)),rep(2,length(bite_time_2)),rep(3,length(bite_time_3)))
                  # Now we sample the first n_hooks*n_hook_sat_level unadjusted
                  all_times <- c(bite_time_1,bite_time_2,bite_time_3)
                  time_ind <- sort.int(all_times, index.return = T, decreasing = F)$ix
                  # for the remaining hooks we sample/thin according to the sat_fun
                  current_sat_level <- hook_sat_level
                  # keep track of the number of hooks currently removed
                  counter <- round(n_hooks*hook_sat_level)
                  for(k2 in (round(n_hooks*hook_sat_level)+1):(length(all_times)))
                  {
                    # create a flag which is true until the time is matched to the correct species group
                    flag <- T
                    if(species_ind[time_ind[k2]]==1)
                    {
                      # Bernoulli trial to determine if the fish is captured
                      # If captured, keep the time, else set to NA
                      time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                    prob=sat_fun(sat_effect = saturation_effect[1],
                                                                 sat_level = current_sat_level,
                                                                 hook_sat_level = hook_sat_level))==1,
                                             time_ind[k2], NA)
                      flag <- F
                    }
                    if(species_ind[time_ind[k2]]==2 & flag)
                    {
                      # Bernoulli trial to determine if the fish is captured
                      time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                    prob=sat_fun(sat_effect = saturation_effect[2],
                                                                 sat_level = current_sat_level,
                                                                 hook_sat_level = hook_sat_level))==1,
                                             time_ind[k2], NA)
                      flag <- F
                    }
                    if(species_ind[time_ind[k2]]==3 & flag)
                    {
                      # Bernoulli trial to determine if the fish is captured
                      time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                    prob=sat_fun(sat_effect = saturation_effect[3],
                                                                 sat_level = current_sat_level,
                                                                 hook_sat_level = hook_sat_level))==1,
                                             time_ind[k2], NA)
                    }
                    if(!is.na(time_ind[k2]))
                    {
                      # update counter to keep track of number of fish caught
                      counter <- counter + 1
                      current_sat_level <- counter/n_hooks
                    }
                    if(counter==n_hooks)
                    {
                      # if all hooks bitten, break the loop
                      time_ind <- time_ind[1:k2]
                      break
                    }
                  }
                  # remove the NA times (keep only the successful captures)
                  time_ind <- time_ind[!is.na(time_ind)]
                  nbite$bites[nbite$species==1 &
                                nbite$station==i2 &
                                nbite$year==j2] <- sum(species_ind[time_ind]==1)
                  nbite$bites[nbite$species==2 &
                                nbite$station==i2 &
                                nbite$year==j2] <- sum(species_ind[time_ind]==2)
                  nbite$bites[nbite$species==3 &
                                nbite$station==i2 &
                                nbite$year==j2] <- sum(species_ind[time_ind]==3)
                }
                
              }
            }
            
            nbite <-
              nbite %>%
              group_by(station, year) %>%
              mutate(prop_sat=sum(bites/n_hooks),
                     composition=bites/(sum(bites)))
            
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
            Results[results_mapper(nsim,i,j,k,l,'censored'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
            Results[results_mapper(nsim,i,j,k,l,'censored'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
            Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
            Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
            
            # fit a naive CPUE model that ignores competition
            dat <- nbite[nbite$species == 3,]
            dat$event_ID <- 1:dim(dat)[1]
            
            # We fit the CPUE method with a parametrization to get at relative abundance
            mod2 <- inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                         data = dat,
                         num.threads = 1, family = 'poisson')
            parameters <- t(sapply(mod2$marginals.fixed, FUN = function(x){
              inla.qmarginal(p=c(0.025, 0.5, 0.975),
                             marginal=inla.tmarginal(fun=exp,marginal=x))}))[-1,]
            colnames(parameters)=c('LCL', 'Median','UCL')
            
            # These measures (Rel_...) are used throughout and are the target
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_Bias'][-1] <-
              (sapply(mod2$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_RMSE'][-1] <-
              (sapply(mod2$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])^2
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_Coverage'][-1] <-
              ifelse(parameters[,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[,3] >= mean_bite[-1,3]/mean_bite[1,3],
                     1,0)
            
            # Next, compute the ICR method. First, scale the catch counts
            dat$bites <- round(dat$bites*comp_factor_fun(1-dat$prop_sat, rep(n_hooks,length(dat$prop_sat))))
            
            # Again - check the model's ability to infer relative abundance (the target)
            mod4 <- inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                         data = dat,
                         num.threads = 1, family = 'poisson')
            
            parameters <- t(sapply(mod4$marginals.fixed, FUN = function(x){
              inla.qmarginal(p=c(0.025, 0.5, 0.975),
                             marginal=inla.tmarginal(fun=exp,marginal=x))}))[-1,]
            colnames(parameters)=c('LCL', 'Median','UCL')
            
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_Bias'][-1] <-
              (sapply(mod4$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_RMSE'][-1] <-
              (sapply(mod4$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])^2
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_Coverage'][-1] <-
              ifelse(parameters[,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[,3] >= mean_bite[-1,3]/mean_bite[1,3],
                     1,0)
            
            # Fit the censored method with hatp*=0.85
            upper_bound <- rep(0, length(nbite[nbite$species==3,]$prop_sat))
            if(cprop < 1)
            {
              # Use the Baranov Catch equation to (starting at 85% saturation) to derive upper bound
              dat <- nbite[nbite$species==3,]
              scale_fac <- rep(0, length(dat$bites))
              scale_fac[dat$prop_sat>cprop] <-
                comp_factor_fun(1-signif((dat[dat$prop_sat>cprop,]$prop_sat-cprop)/(1-cprop),5),
                                rep(round((1-cprop)*n_hooks),sum(dat$prop_sat>cprop)))
              
              upper_bound[dat$prop_sat>cprop] <- round(
                (dat$prop_sat[dat$prop_sat>cprop]-cprop)*n_hooks*
                  scale_fac[dat$prop_sat>cprop])
              
              # How many excess baits could have been removed by the Target species (upper bound)?
              # This will be added to the observed catch count later
              upper_bound <- round(upper_bound)
            }
            
            # Extract species 3 catch counts and fix censorship interval
            dat <- nbite[nbite$species==3,]
            # If low=Inf then INLA treats the data point as uncensored.
            dat$low <- rep(Inf,dim(dat)[1])
            dat$high <- rep(Inf,dim(dat)[1])
            
            # Data points below the max catch count and with p_it > cprop are censored
            dat$low[which(dat$prop_sat >= cprop &
                            0 < upper_bound &
                            dat$bites < quantile(dat$bites,1))] <-
              as.matrix(dat[which(dat$prop_sat >= cprop &
                                    0 < upper_bound &
                                    dat$bites < quantile(dat$bites,1)),
                            c('bites')])[,1]
            
            # For the upper limit, add the ICR-based 'upper_bound' to catch counts
            ## WE HAVE REMOVED THE UPPER BOUND IN THIS EXPERIMENT
            # dat$high[which(dat$prop_sat >= cprop &
            #                  0 < upper_bound &
            #                  dat$bites < quantile(dat$bites,1))] <-
            #   dat$bites[which(dat$prop_sat >= cprop &
            #                     0 < upper_bound &
            #                     dat$bites < quantile(dat$bites,1))] +
            #   upper_bound[which(dat$prop_sat >= cprop &
            #                       0 < upper_bound &
            #                       dat$bites < quantile(dat$bites,1))]
            
            ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))
            
            dat$event_ID <- 1:dim(dat)[1]
            
            # infer the ability to estimate the relative abundance (the target)
            mod6 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                         family="cenpoisson2",verbose=F,
                         data= dat,
                         num.threads = 1,
                         control.mode = list(result=mod2, restart=T))
            parameters <- t(sapply(mod6$marginals.fixed, FUN = function(x){
              inla.qmarginal(p=c(0.025, 0.5, 0.975),
                             marginal=inla.tmarginal(fun=exp,marginal=x))}))[-1,]
            colnames(parameters)=c('LCL', 'Median','UCL')
            
            Results[results_mapper(nsim,i,j,k,l,'censored'),'Rel_Bias'][-1] <-
              (sapply(mod6$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])
            Results[results_mapper(nsim,i,j,k,l,'censored'),'Rel_RMSE'][-1] <-
              (sapply(mod6$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])^2
            Results[results_mapper(nsim,i,j,k,l,'censored'),'Rel_Coverage'][-1] <-
              ifelse(parameters[,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[,3] >= mean_bite[-1,3]/mean_bite[1,3],
                     1,0)
            
            # Repeat for second censored model with hatp* = 0.95
            
            # Need to redefine the upper bound
            upper_bound2 <- rep(0, length(nbite[nbite$species==3,]$prop_sat))
            if(cprop2 < 1)
            {
              # Use the Baranov Catch equation to (starting at 95% saturation) to derive scale factor
              dat <- nbite[nbite$species==3,]
              scale_fac2 <- rep(0, length(dat$bites))
              scale_fac2[dat$prop_sat>cprop2] <-
                comp_factor_fun(1-signif((dat[dat$prop_sat>cprop2,]$prop_sat-cprop2)/(1-cprop2),5),
                                rep(round((1-cprop2)*n_hooks),sum(dat$prop_sat>cprop2)))
              
              upper_bound2[dat$prop_sat>cprop2] <- round(
                (dat$prop_sat[dat$prop_sat>cprop2]-cprop2)*n_hooks*
                  scale_fac[dat$prop_sat>cprop2])
              
              # We will ultimately add this to the catch counts to form the upper bound
              upper_bound2 <- round(upper_bound2)
              
            }
            # Define the censorship intervals
            dat <- nbite[nbite$species==3,]
            dat$low <- rep(Inf,dim(dat)[1])
            dat$high <- rep(Inf,dim(dat)[1])
            
            dat$low[which(dat$prop_sat >= cprop2 &
                            0 < upper_bound2 &
                            dat$bites < quantile(dat$bites,1))] <-
              as.matrix(dat[which(dat$prop_sat >= cprop2 &
                                    0 < upper_bound2 &
                                    dat$bites < quantile(dat$bites,1)),
                            c('bites')])[,1]
            
            ## AGAIN - REMOVE UPPER BOUND
            # dat$high[which(dat$prop_sat >= cprop2 &
            #                  0 < upper_bound2 &
            #                  dat$bites < quantile(dat$bites,1))] <-
            #   dat$bites[which(dat$prop_sat >= cprop2 &
            #                     0 < upper_bound2 &
            #                     dat$bites < quantile(dat$bites,1))] +
            #   upper_bound2[which(dat$prop_sat >= cprop2 &
            #                       0 < upper_bound2 &
            #                       dat$bites < quantile(dat$bites,1))]
            
            ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))
            
            dat$event_ID <- 1:dim(dat)[1]
            
            # Estimate relative abundance (the target)
            mod8 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                         family="cenpoisson2",verbose=F,
                         data= dat,
                         num.threads = 1,
                         control.mode = list(result=mod2, restart=T))
            parameters <- t(sapply(mod8$marginals.fixed, FUN = function(x){
              inla.qmarginal(p=c(0.025, 0.5, 0.975),
                             marginal=inla.tmarginal(fun=exp,marginal=x))}))[-1,]
            colnames(parameters)=c('LCL', 'Median','UCL')
            
            Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Rel_Bias'][-1] <-
              (sapply(mod8$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])
            Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Rel_RMSE'][-1] <-
              (sapply(mod8$marginals.fixed, FUN = function(x){
                inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))})[-1] -
                 mean_bite[-1,3]/mean_bite[1,3])^2
            Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Rel_Coverage'][-1] <-
              ifelse(parameters[,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[,3] >= mean_bite[-1,3]/mean_bite[1,3],
                     1,0)
            
            print(Results[results_mapper(nsim,i,j,k,l,'naive'),])
            print(Results[results_mapper(nsim,i,j,k,l,'adjust'),])
            print(Results[results_mapper(nsim,i,j,k,l,'censored'),])
            print(Results[results_mapper(nsim,i,j,k,l,'censored_95'),])
            
          }
          
        }
      }
    }
  }
  # Save the results
  r_num <- sample(1e6, size = 1)
  saveRDS(Results, paste0('Simulation_Results_Bonus_Corrected_',r_num,'.rds'))
  
}
if(!run_sim)
{
  # Old code used to compile lots of smaller simulation files 
  # speeds up process as can run experiment in batches of 10
  
  # files <- list.files(pattern='Simulation_Results_Bonus')
  # count <- 1
  # for(i in files)
  # {
  #   if(count == 1)
  #   {
  #     Results <- readRDS(i)
  #   }
  #   if(count > 1)
  #   {
  #     Results_tmp <- readRDS(i)
  #     Results_tmp$nsim <- Results_tmp$nsim + (count - 1)*10
  #     Results <- rbind(Results, Results_tmp)
  #   }
  #   count <- count + 1
  # }
  # rm(Results_tmp)
  
  # Can simply read this single file as have compiled them together
  # manually
  Results <- readRDS('Simulation_Results_Bonus_Corrected.rds')
  
  library(inlabru)
  
  # Change the factors to ordered factors to improve plotting
  Results$sat_level <- factor(Results$sat_level, levels=c('low','high'), ordered = T)
  Results$mean_attract <- factor(Results$mean_attract, levels=c('increase','decrease'), ordered = T)
  Results$model <- factor(Results$model, levels=c('naive','adjust','censored','censored_95'), ordered = T)
  
  # Create artificial 'relative abundance' of target and Competitor plots
  rel_abund_dat <- data.frame(expand.grid(
    species=c('Target species','Competitor'),
    sat_level=factor(c('low','high'), levels=c('low','high'), ordered = T),
    mean_attract=factor(c('increase','increase','decrease','decrease'), levels=c('increase','decrease'), ordered = T),
    Year=c(1,2,3,4,5,6)))
  rel_abund_dat$Abundance <- 1
  rel_abund_dat$Abundance[rel_abund_dat$species=='Target species'] <-
    7 - rel_abund_dat$Year[rel_abund_dat$species=='Target species']
  rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                            rel_abund_dat$sat_level=='low'&
                            rel_abund_dat$mean_attract=='increase'] <-
    (c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                                                              rel_abund_dat$sat_level=='low'&
                                                              rel_abund_dat$mean_attract=='increase']]
  rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                            rel_abund_dat$sat_level=='high'&
                            rel_abund_dat$mean_attract=='increase'] <-
    2*(c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                                                                rel_abund_dat$sat_level=='high'&
                                                                rel_abund_dat$mean_attract=='increase']]
  rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                            rel_abund_dat$sat_level=='low'&
                            rel_abund_dat$mean_attract=='decrease'] <-
    (c(280,200,180,160,140,120)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                                                          rel_abund_dat$sat_level=='low'&
                                                          rel_abund_dat$mean_attract=='decrease']]
  rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                            rel_abund_dat$sat_level=='high'&
                            rel_abund_dat$mean_attract=='decrease'] <-
    2*(c(280,200,180,160,140,120)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                                                            rel_abund_dat$sat_level=='high'&
                                                            rel_abund_dat$mean_attract=='decrease']]
  
  rel_abund_plot <-
    ggplot(rel_abund_dat,
           aes(x=Year, y=Abundance, linetype=species) ) +
    geom_rect(data= ~.x[.x$Year==1,],
              aes(x=Year, y=Abundance, linetype=species, fill = mean_attract),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3) +
    geom_line() +
    facet_grid(mean_attract+sat_level~.,
               labeller = labeller(
                 mean_attract = c(
                   increase = 'increasing abundance of \nmain competitor species',
                   decrease = 'decreasing abundance of \nmain competitor species'
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 ))) +
    ylab('')  +
    scale_fill_brewer(palette = 'Pastel2') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank()) +
    guides(fill='none') +
    theme(strip.background = element_blank(),strip.text = element_blank(),
          legend.position = c(0.70,0.97),legend.box.background=element_blank(),
          legend.background=element_blank(),
          axis.title.y = element_blank()) +
    guides(linetype=guide_legend('')) +
    ylim(c(0,6))
  
  # THESE ARE DESIGNED FOR A4 LANDSCAPE or 5.83 x 11.3
  # NOTICE THE HACK IN MULTIPLOT'S LAYOUT ARGUMENT
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1) %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Mean = median(Rel_Bias),
             UCL = quantile(Rel_Bias, prob=0.975),
             LCL = quantile(Rel_Bias, prob=0.025)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0) +
      ggtitle('Bias in relative abundance for each method')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('Bias') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                         values=c('circle','triangle','square','square')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),
             shape=guide_legend(title='Method')),
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1) %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Mean = median(Rel_Bias),
             UCL = median(Rel_Bias)+(2/sqrt(length(Rel_Bias)))*mad(Rel_Bias),
             LCL = median(Rel_Bias)-(2/sqrt(length(Rel_Bias)))*mad(Rel_Bias)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      #geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0) +
      ggtitle('Bias in relative abundance for each method')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('Bias') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                         values=c('circle','triangle','square','square')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),
             shape=guide_legend(title='Method')),
    rel_abund_plot,
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1, !(model %in% c('naive','adjust'))) %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Mean = median(Rel_Bias),
             UCL = median(Rel_Bias)+(2/sqrt(length(Rel_Bias)))*mad(Rel_Bias),
             LCL = median(Rel_Bias)-(2/sqrt(length(Rel_Bias)))*mad(Rel_Bias)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0) +
      ggtitle('Bias in relative abundance for each method')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('Bias') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('Censored 0.85','Censored 0.95'),
                         values=c('circle','triangle')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),          shape=guide_legend(title='Method')),
    rel_abund_plot,
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1) %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Coverage_Mean = mean(Rel_Coverage),
             UCL=mean(Rel_Coverage)+(2/sqrt(length(Rel_Coverage)))*sqrt(mean(Rel_Coverage)*(1-mean(Rel_Coverage))),
             LCL=mean(Rel_Coverage)-(2/sqrt(length(Rel_Coverage)))*sqrt(mean(Rel_Coverage)*(1-mean(Rel_Coverage)))) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Coverage_Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Coverage_Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0.95) +
      ggtitle('Coverage of intervals of relative abundance for each method')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('Coverage') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                         values=c('circle','triangle','square','square')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),          shape=guide_legend(title='Method')),
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1) %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Mean = median(Rel_RMSE),
             UCL = median(Rel_RMSE)+(2/sqrt(length(Rel_RMSE)))*mad(Rel_RMSE),
             LCL = median(Rel_RMSE)-(2/sqrt(length(Rel_RMSE)))*mad(Rel_RMSE)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_RMSE)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0) +
      ggtitle('MSE in relative abundance for each method')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('MSE') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                         values=c('circle','triangle','square','square')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),          shape=guide_legend(title='Method')),
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1, model!='naive') %>%
      group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Mean = median(Rel_RMSE),
             UCL = median(Rel_RMSE)+(2/sqrt(length(Rel_RMSE)))*mad(Rel_RMSE),
             LCL = median(Rel_RMSE)-(2/sqrt(length(Rel_RMSE)))*mad(Rel_RMSE)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_RMSE)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      geom_hline(yintercept=0) +
      ggtitle('MSE in relative abundance for each method') +#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('MSE') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('ICR','Censored 0.85','Censored 0.95')) +
      scale_shape_manual(labels=c('ICR','Censored 0.85','Censored 0.95'),
                         values=c('triangle','square','square')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),          shape=guide_legend(title='Method')),
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
  
  multiplot(
    rel_abund_plot + ggtitle('Simulated abundance'),
    Results %>%
      filter(Station>1) %>%
      group_by(Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
      mutate(Sat85_Mean = mean(Prop_Sat_85),
             Sat100_Mean = mean(Prop_Sat_100)) %>%
      ungroup(Station) %>%
      mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
      group_by(Station) %>%
      ggplot(aes(x=Station, y=Sat85_Mean, colour='85% Saturation')) +
      geom_rect(data= ~.x[.x$random_samp==T,],
                aes(x=Station, y=Sat85_Mean, colour='85% Saturation', fill = mean_attract),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 1) +
      #geom_errorbar(position = position_dodge(width=0.8)) +
      geom_point(position = position_dodge(width=0.8), size=2) +
      geom_point(aes(x=Station, y=Sat100_Mean, colour='100% Saturation'),position = position_dodge(width=0.8), size=2) +
      ylim(c(0,1)) +
      facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
                 labeller = labeller(
                   sat_effects=c(
                     `no saturation` = 'p* = 1',
                     saturation = 'p* = 0.85'
                   ),
                   bite_fun=c(
                     constant = 'Identical Aggressiveness',
                     mixed = 'Differing Aggressiveness'
                   ),
                   mean_attract = c(
                     increase = 'increasing abundance of \nmain competitor species',
                     decrease = 'decreasing abundance of \nmain competitor species '
                   ),
                   sat_level = c(
                     low = 'saturation less "common"',
                     high = 'saturation "common"'
                   )
                 )) +
      ggtitle('Proportion of Fishing Events with Specified Level of Hook Saturation')+#,
      #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
      ylab('Proportion of fishing events') +
      xlab('Year') + guides(fill='none') +
      scale_fill_brewer(palette = 'Pastel2') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            legend.position = 'left', strip.text.y = element_blank()) +
      scale_color_viridis_d(labels=c('100%','>85%')) +
      guides(color=guide_legend(override.aes=list(fill=NA),title='Method'),          shape=guide_legend(title='Method')) +
      labs(colour='Baits Removed'),
    layout = matrix(c(rep(NA,55),rep(1,445),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
  
}  
