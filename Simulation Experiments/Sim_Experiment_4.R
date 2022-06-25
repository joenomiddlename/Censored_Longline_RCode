# Simulation study demonstrating the censored hook competition method
# 3 species with settings the same as experiment 2 but with true p* value equal 0.5
# This code is commented much less as it is the same as experiment 2
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)

nspecies <- 3
bite_funs <- cbind(rep('constant',3),
                   c('exp_decay','mixed', 'constant'))
soak_time <- 5
n_hooks <- 800
# definition of stations and years are opposite to the manuscript
nstation <- 6
nyears <- 30
saturation_level <- c(1,2)
mean_attract <- c('constant', 'linear')
sd_log_bite <- c(0.8, 0.2, 0)
n_sim <- 100
hook_sat_level <- 0.5 # true proportion at which saturation effects begin
cprop=0.85 # assumed proportion at which saturation effects begin
cprop2=0.95 # assumed proportion "" "" second model

# how much does the bite rate of each species (linearly) decrease at 100% saturation
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

# try the hook competition adjusted method
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

for(nsim in 68:n_sim)
{
  print(paste0('iteration ',nsim,' out of ',n_sim))
  for(i in 1:length(saturation_level))
  {
    for(j in 1:length(mean_attract))
    {
      for(k in 1:dim(saturation_effects)[2])
      {
        if(mean_attract[j] == 'constant')
        {
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  rep(100, 6))
        }
        if(mean_attract[j] == 'linear')
        {
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  c(120, 140, 160, 180, 200, 220)-100)
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
              bite_time_1 <- bite_samp(bite_funs[1,l],sum(nbite$attracted[nbite$species==1 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_1)>soak_time)
              {
                bite_time_1[bite_time_1>soak_time] <-
                  bite_samp(bite_funs[1,l],sum(bite_time_1>soak_time))
              }
              # repeat for species 2 from uniform distribution
              bite_time_2 <- bite_samp(bite_funs[2,l],sum(nbite$attracted[nbite$species==2 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_2)>soak_time)
              {
                bite_time_2[bite_time_2>soak_time] <-
                  bite_samp(bite_funs[2,l],sum(bite_time_2>soak_time))
              }

              # repeat for species 3 from uniform distribution
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
                counter <- round(n_hooks*hook_sat_level)
                for(k2 in (round(n_hooks*hook_sat_level)+1):(length(all_times)))
                {
                  flag <- T
                  if(species_ind[time_ind[k2]]==1)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[1],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                    flag <- F
                  }
                  if(species_ind[time_ind[k2]]==2 & flag)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[2],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                    flag <- F
                  }
                  if(species_ind[time_ind[k2]]==3 & flag)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[3],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                  }
                  if(!is.na(time_ind[k2]))
                  {
                    counter <- counter + 1
                    current_sat_level <- counter/n_hooks
                  }
                  if(counter==n_hooks)
                  {
                    time_ind <- time_ind[1:k2]
                    break
                  }
                }
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


          # fit a naive model that ignores competition
          dat <- nbite[nbite$species == 3,]
          dat$event_ID <- 1:dim(dat)[1]
          mod <- inla(bites ~ -1 + factor(station) + f(event_ID, constr=T, model='iid'),
                      data = dat, family = 'poisson')

          parameters <- t(sapply(mod$marginals.fixed, FUN = function(x){
            inla.qmarginal(p=c(0.025, 0.5, 0.975),
                           marginal=inla.tmarginal(fun=exp,marginal=x))}))
          colnames(parameters)=c('LCL', 'Median','UCL')

          Results[results_mapper(nsim,i,j,k,l,'naive'),'Bias'] <-
            (sapply(mod$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])
          Results[results_mapper(nsim,i,j,k,l,'naive'),'RMSE'] <-
            (sapply(mod$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])^2
          Results[results_mapper(nsim,i,j,k,l,'naive'),'Coverage'] <-
            ifelse(parameters[,1] <= mean_bite[,3] & parameters[,3] >= mean_bite[,3],
                   1,0)

          mod2 <- inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                       data = dat, family = 'poisson')

          parameters <- t(sapply(mod2$marginals.fixed, FUN = function(x){
            inla.qmarginal(p=c(0.025, 0.5, 0.975),
                           marginal=inla.tmarginal(fun=exp,marginal=x))}))[-1,]
          colnames(parameters)=c('LCL', 'Median','UCL')

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

          # ICR method
          dat$bites <- round(dat$bites*comp_factor_fun(1-dat$prop_sat, rep(n_hooks,length(dat$prop_sat))))
          mod3 <- inla(bites ~ -1 + factor(station) + f(event_ID, constr=T, model='iid'),
                       data = dat, family = 'poisson')

          parameters <- t(sapply(mod3$marginals.fixed, FUN = function(x){
            inla.qmarginal(p=c(0.025, 0.5, 0.975),
                           marginal=inla.tmarginal(fun=exp,marginal=x))}))
          colnames(parameters)=c('LCL', 'Median','UCL')

          Results[results_mapper(nsim,i,j,k,l,'adjust'),'Bias'] <-
            (sapply(mod3$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])
          Results[results_mapper(nsim,i,j,k,l,'adjust'),'RMSE'] <-
            (sapply(mod3$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])^2
          Results[results_mapper(nsim,i,j,k,l,'adjust'),'Coverage'] <-
            ifelse(parameters[,1] <= mean_bite[,3] & parameters[,3] >= mean_bite[,3],
                   1,0)

          mod4 <- inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                       data = dat, family = 'poisson')

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

          # Censored methods

          # define upper bound
          upper_bound <- rep(0, length(nbite[nbite$species==3,]$prop_sat))
          species_1_uncensoredprop <- 0.5
          if(cprop < 1)
          {
            # Use the Baranov Catch equation to (starting at 85% saturation) to derive scale factor
            dat <- nbite[nbite$species==3,]
            scale_fac <- rep(0, length(dat$bites))
            scale_fac[dat$prop_sat>cprop] <-
              comp_factor_fun(1-signif((dat[dat$prop_sat>cprop,]$prop_sat-cprop)/(1-cprop),5),
                              rep(round((1-cprop)*n_hooks),sum(dat$prop_sat>cprop)))

            upper_bound[dat$prop_sat>cprop] <- round(
              (dat$prop_sat[dat$prop_sat>cprop]-cprop)*n_hooks*
                scale_fac[dat$prop_sat>cprop])

            # How many excess hooks could have gone to the Target species?
            upper_bound <- round(upper_bound)
          }
          dat <- nbite[nbite$species==3,]
          dat$low <- rep(Inf,dim(dat)[1])
          dat$high <- rep(Inf,dim(dat)[1])

          dat$low[which(dat$prop_sat >= cprop &
                          0 < upper_bound &
                          dat$bites < quantile(dat$bites,1))] <-
            as.matrix(dat[which(dat$prop_sat >= cprop &
                                  0 < upper_bound &
                                  dat$bites < quantile(dat$bites,1)),
                          c('bites')])[,1]

          dat$high[which(dat$prop_sat >= cprop &
                           0 < upper_bound &
                           dat$bites < quantile(dat$bites,1))] <-
            dat$bites[which(dat$prop_sat >= cprop &
                              0 < upper_bound &
                              dat$bites < quantile(dat$bites,1))] +
            upper_bound[which(dat$prop_sat >= cprop &
                                0 < upper_bound &
                                dat$bites < quantile(dat$bites,1))]


          ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

          dat$event_ID <- 1:dim(dat)[1]
          # Attempt to estimate absolute abundance (futile)
          mod5 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ -1 + factor(station) + f(event_ID, constr=T, model='iid'),
                       family="cenpoisson2",verbose=F,
                       data= dat)

          parameters <- t(sapply(mod5$marginals.fixed, FUN = function(x){
            inla.qmarginal(p=c(0.025, 0.5, 0.975),
                           marginal=inla.tmarginal(fun=exp,marginal=x))}))
          colnames(parameters)=c('LCL', 'Median','UCL')

          Results[results_mapper(nsim,i,j,k,l,'censored'),'Bias'] <-
            (sapply(mod5$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])
          Results[results_mapper(nsim,i,j,k,l,'censored'),'RMSE'] <-
            (sapply(mod5$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])^2
          Results[results_mapper(nsim,i,j,k,l,'censored'),'Coverage'] <-
            ifelse(parameters[,1] <= mean_bite[,3] & parameters[,3] >= mean_bite[,3],
                   1,0)
          # Attempt to estimate relative abundance (target)
          mod6 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                       family="cenpoisson2",verbose=F,
                       data= dat)

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

          # Repeat for different hatp* value

          # censorship interval
          upper_bound2 <- rep(0, length(nbite[nbite$species==3,]$prop_sat))
          if(cprop2 < 1)
          {
            # Use the Baranov Catch equation to derive upper bound
            dat <- nbite[nbite$species==3,]
            scale_fac2 <- rep(0, length(dat$bites))
            scale_fac2[dat$prop_sat>cprop2] <-
              comp_factor_fun(1-signif((dat[dat$prop_sat>cprop2,]$prop_sat-cprop2)/(1-cprop2),5),
                              rep(round((1-cprop2)*n_hooks),sum(dat$prop_sat>cprop2)))

            upper_bound2[dat$prop_sat>cprop2] <- round(
              (dat$prop_sat[dat$prop_sat>cprop2]-cprop2)*n_hooks*
                scale_fac[dat$prop_sat>cprop2])

            upper_bound2 <- round(upper_bound2)
          }
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

          dat$high[which(dat$prop_sat >= cprop2 &
                           0 < upper_bound2 &
                           dat$bites < quantile(dat$bites,1))] <-
            dat$bites[which(dat$prop_sat >= cprop2 &
                              0 < upper_bound2 &
                              dat$bites < quantile(dat$bites,1))] +
            upper_bound2[which(dat$prop_sat >= cprop2 &
                                0 < upper_bound2 &
                                dat$bites < quantile(dat$bites,1))]


          ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

          dat$event_ID <- 1:dim(dat)[1]

          # We try to estimate the absolute abundance (futile)
          mod7 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ -1 + factor(station) + f(event_ID, constr=T, model='iid'),
                       family="cenpoisson2",verbose=F,
                       data= dat)

          parameters <- t(sapply(mod7$marginals.fixed, FUN = function(x){
            inla.qmarginal(p=c(0.025, 0.5, 0.975),
                           marginal=inla.tmarginal(fun=exp,marginal=x))}))
          colnames(parameters)=c('LCL', 'Median','UCL')

          Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Bias'] <-
            (sapply(mod7$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])
          Results[results_mapper(nsim,i,j,k,l,'censored_95'),'RMSE'] <-
            (sapply(mod7$marginals.fixed, FUN = function(x){
              inla.emarginal(fun=function(x){return(x)},marginal=inla.tmarginal(fun=exp,marginal=x))}) -
               mean_bite[,3])^2
          Results[results_mapper(nsim,i,j,k,l,'censored_95'),'Coverage'] <-
            ifelse(parameters[,1] <= mean_bite[,3] & parameters[,3] >= mean_bite[,3],
                   1,0)

          # We try to estimate the relative abundance (target)
          mod8 <- inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                       family="cenpoisson2",verbose=F,
                       data= dat)

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
# Save results
#saveRDS(Results, 'Simulation_Results_Corrected_cen50.rds')
Results <- readRDS('Simulation_Results_Corrected_cen50.rds')

library(inlabru)

# Change the factors to ordered factors to improve plotting
Results$sat_level <- factor(Results$sat_level, levels=c('low','high'), ordered = T)
Results$mean_attract <- factor(Results$mean_attract, levels=c('constant','linear'), ordered = T)
Results$model <- factor(Results$model, levels=c('naive','adjust','censored','censored_95'), ordered = T)

# Create artificial 'relative abundance' of target and Competitor plots
rel_abund_dat <- data.frame(expand.grid(
  species=c('Target species','Competitor'),
  sat_level=factor(c('low','high'), levels=c('low','high'), ordered = T),
  mean_attract=factor(c('constant','constant','linear','linear')),
  Year=c(1,2,3,4,5,6)))
rel_abund_dat$Abundance <- 1
rel_abund_dat$Abundance[rel_abund_dat$species=='Target species'&
                          rel_abund_dat$mean_attract=='linear'] <-
  rel_abund_dat$Year[rel_abund_dat$species=='Target species'&
                       rel_abund_dat$mean_attract=='linear']
rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                          rel_abund_dat$sat_level=='low'] <-
  (c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                       rel_abund_dat$sat_level=='low']]
rel_abund_dat$Abundance[rel_abund_dat$species=='Competitor'&
                          rel_abund_dat$sat_level=='high'] <-
  2*(c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='Competitor'&
                       rel_abund_dat$sat_level=='high']]

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
                 constant = 'constant Target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               ))) +
  ylab('')  +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  guides(fill='none') +
  theme(strip.background = element_blank(),strip.text = element_blank(),
        legend.position = c(0.43,0.94),legend.box.background=element_blank(),
        legend.background=element_blank(),
        axis.title.y = element_blank()) +
  guides(linetype=guide_legend('')) +
  ylim(c(0,6))

# THESE ARE DESIGNED FOR 5.83x11.3
# NOTICE THE HACK IN MULTIPLOT'S LAYOUT ARGUMENT
multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
Results %>%
  filter(Station>1, sat_effects!='no saturation') %>%
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
  facet_grid(mean_attract + sat_level ~ bite_fun , scales = 'free_y',
             labeller = labeller(

               bite_fun=c(
                 constant = 'Identical Aggressiveness',
                 mixed = 'Differing Aggressiveness'
               ),
               mean_attract = c(
                 constant = 'constant Target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               )
             )) +
  geom_hline(yintercept=0) +
  ggtitle('Bias in relative abundance for each method')+#,
          #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
  ylab('Bias') +
  xlab('Year') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        legend.position = 'left', strip.text.y = element_blank()) +
  scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
  scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                     values=c('circle','triangle','square','square')) +
  guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
         shape=guide_legend(title='Method')),
rel_abund_plot,
layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
Results %>%
  filter(Station>1, sat_effects!='no saturation') %>%
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
                 saturation = 'p* = 0.5'
               ),
               bite_fun=c(
                 constant = 'Identical Aggressiveness',
                 mixed = 'Differing Aggressiveness'
               ),
               mean_attract = c(
                 constant = 'constant Target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               )
             )) +
  geom_hline(yintercept=0) +
  ggtitle('Bias in relative abundance for each method')+#,
          #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
  ylab('Bias') +
  xlab('Year') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        legend.position = 'left', strip.text.y = element_blank()) +
  scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
  scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                     values=c('circle','triangle','square','square')) +
  guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
         shape=guide_legend(title='Method')),
layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
Results %>%
  filter(Station>1, !(model %in% c('naive','adjust')), sat_effects!='no saturation') %>%
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
                 saturation = 'p* = 0.5'
               ),
               bite_fun=c(
                 constant = 'Identical Aggressiveness',
                 mixed = 'Differing Aggressiveness'
               ),
               mean_attract = c(
                 constant = 'constant Target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               )
             )) +
  geom_hline(yintercept=0) +
  ggtitle('Bias in relative abundance for each method')+#,
          #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
  ylab('Bias') +
  xlab('Year') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        legend.position = 'left', strip.text.y = element_blank()) +
  scale_color_viridis_d(labels=c('Censored 85','Censored 95')) +
  scale_shape_manual(labels=c('Censored 85','Censored 95'),
                     values=c('circle','triangle')) +
  guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
         shape=guide_legend(title='Method')),
layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
Results %>%
  filter(Station>1, sat_effects!='no saturation') %>%
  group_by(model, Station, bite_fun,  sat_level, sat_effects, mean_attract) %>%
  mutate(Coverage_Mean = mean(Rel_Coverage),
         UCL=mean(Rel_Coverage)+(2/sqrt(100))*sqrt(mean(Rel_Coverage)*(1-mean(Rel_Coverage))),
         LCL=mean(Rel_Coverage)-(2/sqrt(100))*sqrt(mean(Rel_Coverage)*(1-mean(Rel_Coverage)))) %>%
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
                 saturation = 'p* = 0.5'
               ),
               bite_fun=c(
                 constant = 'Identical Aggressiveness',
                 mixed = 'Differing Aggressiveness'
               ),
               mean_attract = c(
                 constant = 'constant Target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               )
             )) +
  geom_hline(yintercept=0.95) +
  ggtitle('Coverage of intervals of relative abundance for each method')+#,
          #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
  ylab('Coverage') +
  xlab('Year') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        legend.position = 'left', strip.text.y = element_blank()) +
  scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
  scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                     values=c('circle','triangle','square','square')) +
  guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
         shape=guide_legend(title='Method')),
layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(Station>1, sat_effects!='no saturation') %>%
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
    #geom_errorbar(position = position_dodge(width=0.8)) +
    geom_point(position = position_dodge(width=0.8), size=2) +
    facet_grid(mean_attract + sat_level ~ bite_fun + sat_effects , scales = 'free_y',
               labeller = labeller(
                 sat_effects=c(
                   `no saturation` = 'p* = 1',
                   saturation = 'p* = 0.5'
                 ),
                 bite_fun=c(
                   constant = 'Identical Aggressiveness',
                   mixed = 'Differing Aggressiveness'
                 ),
                 mean_attract = c(
                   constant = 'constant Target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0) +
    ggtitle('MSE in relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title='Method')),
  layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(Station>1, sat_effects!='no saturation', model!='naive') %>%
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
                   saturation = 'p* = 0.5'
                 ),
                 bite_fun=c(
                   constant = 'Identical Aggressiveness',
                   mixed = 'Differing Aggressiveness'
                 ),
                 mean_attract = c(
                   constant = 'constant Target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0) +
    ggtitle('MSE in relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of Target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 50% of baits removed.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 0.85','Censored 0.95'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title='Method')),
  layout = matrix(c(rep(NA,35),rep(1,465),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))


Results %>%
  group_by(Station, mean_attract,  sat_level) %>%
  #mutate(Mean = mean(Prop_Sat),
  #       UCL = quantile(Prop_Sat, prob=0.975),
  #       LCL = quantile(Prop_Sat, prob=0.025)) %>%
  ggplot(aes(x=Station, y=Prop_Sat_85))+#, ymin=LCL, ymax=UCL)) +
  #geom_errorbar(position = position_dodge(width=0.3)) +
  geom_point(position = position_dodge(width=0.3)) +
  facet_grid(mean_attract ~ sat_level)

Results %>%
  group_by(Station, mean_attract,  sat_level) %>%
  #mutate(Mean = mean(Prop_Sat),
  #       UCL = quantile(Prop_Sat, prob=0.975),
  #       LCL = quantile(Prop_Sat, prob=0.025)) %>%
  ggplot(aes(x=Station, y=Prop_Sat_100))+#, ymin=LCL, ymax=UCL)) +
  #geom_errorbar(position = position_dodge(width=0.3)) +
  geom_point(position = position_dodge(width=0.3)) +
  facet_grid(mean_attract ~ sat_level)

Results %>%
  filter(Station>1, !(model %in% c('naive','adjust') ) )%>%
  group_by(sat_effects, bite_fun, sat_level, mean_attract) %>%
  mutate(Mean = mean(Prop_Sat_85),
         Mean2 = mean(Prop_Sat_100)) %>%
  ggplot() +
  geom_hline(aes(yintercept=Mean), linetype='solid', colour='black') +
  geom_hline(aes(yintercept=Mean2), linetype='dotted', colour='red') +
  facet_grid(sat_level + mean_attract ~ bite_fun + sat_effects) +
  ylab('Convergence Proportion') +
  ggtitle('Proportion of 100 Simulations That Converged for each method')+#,
          #subtitle = 'Rows are degree of saturation and trend in relative abundance,\nColumns are correlation of Target species\' abundance with saturation events \nSolid black line indicates proportion of converged simulations with 85% bait removal \nDotted red line indicates proportion of converged simulations with 100% bait removal') +
  xlab('Year') +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(colour='none')
