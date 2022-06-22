# Simulation study demonstrating the censored hook competition method
# 1 species who self-censors itself due to abundance exceeding n_hooks during some fishing events
# Simulates spiny dogfish scenario where group sizes of many hundreds or thousands occur
# See https://www.dfo-mpo.gc.ca/species-especes/profiles-profils/spiny-dogfish-aiguillat-commun-atl-eng.html
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)

nspecies <- 1
bite_funs <- 'constant'
soak_time <- 5
n_hooks <- 800
# Note that throughout the sim experiments, the labels of stations and years are swapped.
nstation <- 6
nyears <- 30
# Is the abundance of the target species stationary or increasing?
mean_attract <- c('constant', 'linear')
saturation_level <- 1
# SD of the lognormal distribution
sd_log_bite <- c(0.8)
n_sim <- 100
hook_sat_level <- 0.85 # true proportion at which saturation effects begin
cprop=0.85 # assumed proportion at which saturation effects begin
cprop2=0.95 # assumed proportion "" "" second model

# how much does the bite rate of each species (linearly) decrease at 100% saturation
# assume linear decrease from 85% saturation onwards
saturation_effects <- cbind(c(0),
                           c(0.8))
# create saturation function
sat_fun <- function(sat_effect, sat_level, hook_sat_level=0.85)
{
  val <- rep(1, length(sat_level))

  val[which(sat_level>hook_sat_level)] <-
    (1 - sat_effect*((sat_level[which(sat_level>hook_sat_level)]-hook_sat_level)/(1-hook_sat_level)))
  return(val)
}

# the hook competition adjustment function
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

# Create a df for storing the results
Results <- data.frame(
  nsim = rep(1:n_sim, each=nstation*1*2*2*1*4),
  sat_level = rep(rep(c('high'), times=n_sim),each=2*2*1*4*nstation),
  mean_attract = rep(rep(mean_attract, times=n_sim*1), each=2*1*4*nstation),
  sat_effects = rep(rep(c('no saturation','saturation'),times=n_sim*1*2), each=1*4*nstation),
  bite_fun = rep(rep(c('constant'), times=n_sim*1*2*2), each=4*nstation),
  model=rep(rep(c('naive','adjust','censored_cprop1','censored_cprop1_85'), times=n_sim*1*2*2), each=nstation),
  Bias=rep(0, times=n_sim*1*2*2*1*4*nstation),
  RMSE=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Coverage=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Rel_Bias=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Rel_RMSE=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Rel_Coverage=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Station=rep(1:nstation, times=n_sim*1*2*2*1*4),
  Prop_Sat_85=rep(0, times=n_sim*1*2*2*1*4*nstation),
  Prop_Sat_100=rep(0, times=n_sim*1*2*2*1*4*nstation))

# A function for mapping each sim iteration to the correct df row
results_mapper <- function(n,i,j,k,l,mod)
{
  return(
    which(Results$nsim == n & Results$sat_level == c('high')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$sat_effects == c('no saturation','saturation')[k] &
            Results$bite_fun == c('constant')[l] &
            Results$model == mod)
         )
}

for(nsim in 1:n_sim)
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
            cbind(rep(800,6)*saturation_level[i]/sd_log_bite[1])
          mean_bite_gen =
            cbind(rep(800,6)*saturation_level[i])
        }
        if(mean_attract[j] == 'linear')
        {
          mean_bite =
            cbind(c(300, 400, 500, 600, 700, 800)*saturation_level[i]/sd_log_bite[1])
          mean_bite_gen =
            cbind(c(300, 400, 500, 600, 700, 800)*saturation_level[i])
        }
          for(l in 1:length(bite_funs))
          {
          saturation_effect <- saturation_effects[,k]
          # sample the number of each species that WOULD bite at each station for each year if hooks were available
          nbite <- data.frame(bites = rep(0, times=nspecies*nstation*nyears),
                              attracted = rep(0, times=nspecies*nstation*nyears),
                              species=rep(1, each=nstation*nyears),
                              station=rep(1:nstation, times=nspecies*nyears),
                              year=rep(rep(1:nyears, each=nstation), times=nspecies))
          nbite$attracted <- rpois(dim(nbite)[1],
                                   lambda = as.numeric(mean_bite[cbind(nbite$station,nbite$species)])*
                                     exp(rnorm(dim(nbite)[1], mean = 0, sd=sd_log_bite[nbite$species])))

          for(i2 in 1:nstation)
          {
            for(j2 in 1:nyears)
            {
              bite_time_1 <- bite_samp(bite_funs[l],sum(nbite$attracted[nbite$species==1 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_1)>soak_time)
              {
                bite_time_1[bite_time_1>soak_time] <-
                  bite_samp(bite_funs[l],sum(bite_time_1>soak_time))
              }

              # Now we sample the first n_hooks*n_hook_sat_level unadjusted
              if((length(bite_time_1)) <= n_hooks*hook_sat_level)
              {
                nbite$bites[nbite$species==1 &
                              nbite$station==i2 &
                              nbite$year==j2] <- length(bite_time_1)
              }
              if((length(bite_time_1)) > n_hooks*hook_sat_level)
              {
                species_ind <- c(rep(1, length(bite_time_1)))
                # Now we sample the first n_hooks*n_hook_sat_level unadjusted
                all_times <- c(bite_time_1)
                # sort the bite times in increasing order
                time_ind <- sort.int(all_times, index.return = T, decreasing = F)$ix
                # for the remaining hooks we sample/thin according to the sat_fun
                current_sat_level <- hook_sat_level
                # keep track of how many hooks have been taken
                counter <- round(n_hooks*hook_sat_level)
                # loop through the remaining hooks (the ones for which capture prob < 1).
                for(k2 in (round(n_hooks*hook_sat_level)+1):(length(all_times)))
                {
                  if(species_ind[time_ind[k2]]==1) # almost always satisfied
                  {
                    # if bernoulli trial success -> keep bite time, else replace with NA
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect,
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                  }
                  if(!is.na(time_ind[k2]))
                  {
                    # if success, then increment counter and sat level
                    counter <- counter + 1
                    current_sat_level <- counter/n_hooks
                  }
                  if(counter==n_hooks)
                  {
                    # if we have caught n_hooks fish then break the loop
                    # extract the k2 bite times (and NAs)
                    time_ind <- time_ind[1:k2]
                    break
                  }
                }
                # filter out the NAs
                time_ind <- time_ind[!is.na(time_ind)]
                # Add the catch count
                # Note - the species==1 is always true at the moment
                # species variable added to the dataframe to allow easy
                # future modification to record multiple species' catch counts
                nbite$bites[nbite$species==1 &
                              nbite$station==i2 &
                              nbite$year==j2] <- sum(species_ind[time_ind]==1)
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
          Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
          Results[results_mapper(nsim,i,j,k,l,'naive'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
          Results[results_mapper(nsim,i,j,k,l,'adjust'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
          Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
          Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
          Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))

          # fit a naive CPUE model that ignores competition
          dat <- nbite[nbite$species == 1,]
          dat$event_ID <- 1:dim(dat)[1]

          mod2 <-
            tryCatch(
              {
                inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                     data = dat, family = 'poisson',verbose=F, control.fixed = list(prec.intercept=1e-1),
                     control.compute = list(config=T))
              },
              error=function(cond)
              {
                return(NULL)
              },
              finally={

              }
            )
          # Has the model converged?
          if(sum(is.na(mod2$summary.fixed$mean))>0 | sum(mod2$misc$mode.status!=0)>0 )
          {
            mod2 <- NULL
          }
          # If yes, store results and add convergence flag
          if(!is.null(mod2))
          {
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Converge'] <- 1

            # Sample parameters from approx posterior distribution and compute summary stats
            parameters <- inla.posterior.sample(1000,mod2,
                                                selection = list(`(Intercept)` = 0,
                                                                 `factor(station)2` = 0,`factor(station)3` = 0,
                                                                 `factor(station)4` = 0,`factor(station)5` = 0,
                                                                 `factor(station)6` = 0))
            parameters <- inla.posterior.sample.eval(fun=function(...){
              c(exp(`(Intercept)`), exp(`factor(station)2`),
                exp(`factor(station)3`), exp(`factor(station)4`),
                exp(`factor(station)5`), exp(`factor(station)6`),
                exp(`(Intercept)`+`factor(station)2`), exp(`(Intercept)`+`factor(station)3`),
                exp(`(Intercept)`+`factor(station)4`), exp(`(Intercept)`+`factor(station)5`),
                exp(`(Intercept)`+`factor(station)6`))}, parameters)
            parameters <- t(apply(parameters, 1, FUN=function(x){
              return(quantile(x,probs=c(0.025,0.5,0.975)))
            }))
            colnames(parameters)=c('LCL', 'Median','UCL')

            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_Bias'][-1] <-
              parameters[2:6,2] -
              mean_bite[-1,1]/mean_bite[1,1]

            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_RMSE'][-1] <-
              (parameters[2:6,2] -
                 mean_bite[-1,1]/mean_bite[1,1])^2

            Results[results_mapper(nsim,i,j,k,l,'naive'),'Rel_Coverage'][-1] <-
              ifelse(parameters[2:6,1] <= mean_bite[-1,1]/mean_bite[1,1] & parameters[2:6,3] >= mean_bite[-1,1]/mean_bite[1,1],
                     1,0)
            # Compare Medians! Median value equals mean_bite_gen by definition of conditional medians
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Bias'] <-
              parameters[c(1,7:11),2] - mean_bite_gen[,1]
            Results[results_mapper(nsim,i,j,k,l,'naive'),'RMSE'] <-
              (parameters[c(1,7:11),2] - mean_bite_gen[,1])^2
            Results[results_mapper(nsim,i,j,k,l,'naive'),'Coverage'] <-
              ifelse(parameters[c(1,7:11),1] <= mean_bite[,1] & parameters[c(1,7:11),3] >= mean_bite[,1],
                     1,0)
          }

          # Next fit the ICR method - first scale the observations by scale factor
          dat$bites <- round(dat$bites*comp_factor_fun(1-dat$prop_sat, rep(n_hooks,length(dat$prop_sat))))

          mod4 <-
            tryCatch(
              {
                inla(bites ~ factor(station) + f(event_ID, constr=T, model='iid'),
                     data = dat, family = 'poisson',verbose=F,
                     control.fixed = list(prec.intercept=1e-1),
                     control.compute = list(config=T))
              },
              error=function(cond)
              {
                return(NULL)
              },
              finally={

              }
            )
          if(sum(is.na(mod4$summary.fixed$mean))>0 | sum(mod4$misc$mode.status!=0)>0 )
          {
            mod4 <- NULL
          }

          if(!is.null(mod4))
          {
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Converge'] <- 1

            parameters <- inla.posterior.sample(1000,mod4,
                                                selection = list(`(Intercept)` = 0,
                                                                 `factor(station)2` = 0,`factor(station)3` = 0,
                                                                 `factor(station)4` = 0,`factor(station)5` = 0,
                                                                 `factor(station)6` = 0))
            parameters <- inla.posterior.sample.eval(fun=function(...){
              c(exp(`(Intercept)`), exp(`factor(station)2`),
                exp(`factor(station)3`), exp(`factor(station)4`),
                exp(`factor(station)5`), exp(`factor(station)6`),
                exp(`(Intercept)`+`factor(station)2`), exp(`(Intercept)`+`factor(station)3`),
                exp(`(Intercept)`+`factor(station)4`), exp(`(Intercept)`+`factor(station)5`),
                exp(`(Intercept)`+`factor(station)6`))}, parameters)
            parameters <- t(apply(parameters, 1, FUN=function(x){
              return(quantile(x,probs=c(0.025,0.5,0.975)))
            }))
            colnames(parameters)=c('LCL', 'Median','UCL')

            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_Bias'][-1] <-
              parameters[2:6,2] -
              mean_bite[-1,1]/mean_bite[1,1]

            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_RMSE'][-1] <-
              (parameters[2:6,2] -
                 mean_bite[-1,1]/mean_bite[1,1])^2

            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Rel_Coverage'][-1] <-
              ifelse(parameters[2:6,1] <= mean_bite[-1,1]/mean_bite[1,1] & parameters[2:6,3] >= mean_bite[-1,1]/mean_bite[1,1],
                     1,0)
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Bias'] <-
              parameters[c(1,7:11),2] - mean_bite_gen[,1]
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'RMSE'] <-
              (parameters[c(1,7:11),2] - mean_bite_gen[,1])^2
            Results[results_mapper(nsim,i,j,k,l,'adjust'),'Coverage'] <-
              ifelse(parameters[c(1,7:11),1] <= mean_bite[,1] & parameters[c(1,7:11),3] >= mean_bite[,1],
                     1,0)
          }

          # Censored method p*=1 : censorship interval
          dat <- nbite[nbite$species==1,]

          # If 'low' is Inf, then INLA treats data as uncensored
          dat$low <- rep(Inf,dim(dat)[1])
          dat$high <- rep(Inf,dim(dat)[1])

          # Treat all data with p_it = 1 as censored
          dat$low[which(dat$prop_sat == 1)] <-
            as.matrix(dat[which(dat$prop_sat == 1),
                          c('bites')])[,1]

          ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

          dat$event_ID <- 1:dim(dat)[1]

          mod6 <-
            tryCatch(
              {
                inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                     family="cenpoisson2",verbose=F,
                     data= dat, control.fixed = list(prec.intercept=1e-1),
                     control.compute = list(config=T),
                     control.mode = list(result=mod2, restart=T))
              },
              error=function(cond)
              {
                return(NULL)
              },
              finally={

              }
            )
          if(sum(is.na(mod6$summary.fixed$mean))>0 | sum(mod6$misc$mode.status!=0)>0 )
          {
            mod6 <- NULL
          }

          if(!is.null(mod6))
          {
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Converge'] <- 1

            parameters <- inla.posterior.sample(1000,mod6,
                                                selection = list(`(Intercept)` = 0,
                                                                 `factor(station)2` = 0,`factor(station)3` = 0,
                                                                 `factor(station)4` = 0,`factor(station)5` = 0,
                                                                 `factor(station)6` = 0))
            parameters <- inla.posterior.sample.eval(fun=function(...){
              c(exp(`(Intercept)`), exp(`factor(station)2`),
                exp(`factor(station)3`), exp(`factor(station)4`),
                exp(`factor(station)5`), exp(`factor(station)6`),
                exp(`(Intercept)`+`factor(station)2`), exp(`(Intercept)`+`factor(station)3`),
                exp(`(Intercept)`+`factor(station)4`), exp(`(Intercept)`+`factor(station)5`),
                exp(`(Intercept)`+`factor(station)6`))}, parameters)
            parameters <- t(apply(parameters, 1, FUN=function(x){
              return(quantile(x,probs=c(0.025,0.5,0.975)))
            }))
            colnames(parameters)=c('LCL', 'Median','UCL')

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Rel_Bias'][-1] <-
              parameters[2:6,2] -
              mean_bite[-1,1]/mean_bite[1,1]

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Rel_RMSE'][-1] <-
              (parameters[2:6,2] -
                 mean_bite[-1,1]/mean_bite[1,1])^2

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Rel_Coverage'][-1] <-
              ifelse(parameters[2:6,1] <= mean_bite[-1,1]/mean_bite[1,1] & parameters[2:6,3] >= mean_bite[-1,1]/mean_bite[1,1],
                     1,0)
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Bias'] <-
              parameters[c(1,7:11),2] - mean_bite_gen[,1]
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'RMSE'] <-
              (parameters[c(1,7:11),2] - mean_bite_gen[,1])^2
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),'Coverage'] <-
              ifelse(parameters[c(1,7:11),1] <= mean_bite[,1] & parameters[c(1,7:11),3] >= mean_bite[,1],
                     1,0)
          }

          # Repeat for other censored models with differing upper observations

          dat <- nbite[nbite$species==1,]
          dat$low <- rep(Inf,dim(dat)[1])
          dat$high <- rep(Inf,dim(dat)[1])

          # Compute the annual 85th catch count percentile
          quant_regions <-
            as.numeric(by(dat$bites,
                          dat$station,
                          FUN = function(x){quantile(x,0.85, na.rm=T)}, simplify = T))

          quant_regions2 <- rep(0, nstation)

          # Only treat the catch counts below the quantile as censored
          dat$low[which(dat$prop_sat == 1 &
                        dat$bites <= quant_regions[dat$station])] <-
            as.matrix(dat[which(dat$prop_sat == 1 &
                                dat$bites <= quant_regions[dat$station]),
                          c('bites')])[,1]

          ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

          dat$event_ID <- 1:dim(dat)[1]

          mod8 <-
            tryCatch(
              {
                inla(formula = inla.mdata(cbind(bites,low,high)) ~ factor(station) + f(event_ID, constr=T, model='iid'),
                     family="cenpoisson2",verbose=F,
                     data= dat, control.fixed = list(prec.intercept=1e-1),
                     control.compute = list(config=T),
                     control.mode = list(result=mod2, restart=T))
              },
              error=function(cond)
              {
                return(NULL)
              },
              finally={

              }
            )
          if(sum(is.na(mod8$summary.fixed$mean))>0 | sum(mod8$misc$mode.status!=0)>0 )
          {
            mod8 <- NULL
          }

          if(!is.null(mod8))
          {
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Converge'] <- 1

            parameters <- inla.posterior.sample(1000,mod8,
                                                selection = list(`(Intercept)` = 0,
                                                                 `factor(station)2` = 0,`factor(station)3` = 0,
                                                                 `factor(station)4` = 0,`factor(station)5` = 0,
                                                                 `factor(station)6` = 0))
            parameters <- inla.posterior.sample.eval(fun=function(...){
              c(exp(`(Intercept)`), exp(`factor(station)2`),
                exp(`factor(station)3`), exp(`factor(station)4`),
                exp(`factor(station)5`), exp(`factor(station)6`),
                exp(`(Intercept)`+`factor(station)2`), exp(`(Intercept)`+`factor(station)3`),
                exp(`(Intercept)`+`factor(station)4`), exp(`(Intercept)`+`factor(station)5`),
                exp(`(Intercept)`+`factor(station)6`))}, parameters)
            parameters <- t(apply(parameters, 1, FUN=function(x){
              return(quantile(x,probs=c(0.025,0.5,0.975)))
            }))
            colnames(parameters)=c('LCL', 'Median','UCL')

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Rel_Bias'][-1] <-
              parameters[2:6,2] -
              mean_bite[-1,1]/mean_bite[1,1]

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Rel_RMSE'][-1] <-
              (parameters[2:6,2] -
                 mean_bite[-1,1]/mean_bite[1,1])^2

            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Rel_Coverage'][-1] <-
              ifelse(parameters[2:6,1] <= mean_bite[-1,1]/mean_bite[1,1] & parameters[2:6,3] >= mean_bite[-1,1]/mean_bite[1,1],
                     1,0)
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Bias'] <-
              parameters[c(1,7:11),2] - mean_bite_gen[,1]
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'RMSE'] <-
              (parameters[c(1,7:11),2] - mean_bite_gen[,1])^2
            Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),'Coverage'] <-
              ifelse(parameters[c(1,7:11),1] <= mean_bite[,1] & parameters[c(1,7:11),3] >= mean_bite[,1],
                     1,0)
          }

          print(Results[results_mapper(nsim,i,j,k,l,'naive'),])
          print(Results[results_mapper(nsim,i,j,k,l,'adjust'),])
          print(Results[results_mapper(nsim,i,j,k,l,'censored_cprop1'),])
          print(Results[results_mapper(nsim,i,j,k,l,'censored_cprop1_85'),])

        }

      }
    }
  }
}
# Save the results
#saveRDS(Results, 'Simulation_Results_Corrected_Self_Censoring.rds')
Results <- readRDS('Simulation_Results_Corrected_Self_Censoring.rds')

library(inlabru)

# Change the factors to ordered factors to improve plotting
Results$sat_level <- factor(Results$sat_level, levels=c('high'), ordered = T)
Results$mean_attract <- factor(Results$mean_attract, levels=c('constant', 'linear'), ordered = T)
Results$model <- factor(Results$model, levels=c('naive','adjust','censored_cprop1','censored_cprop1_85'), ordered = T)

# Create artificial 'relative abundance' of target and aggressive species plots
rel_abund_dat <- data.frame(expand.grid(
  species=c('target species'),
  sat_level=factor(c('high'), levels=c('high'), ordered = T),
  mean_attract=factor(c('constant', 'linear')),
  Year=c(1,2,3,4,5,6)))
rel_abund_dat$Abundance <- 800
rel_abund_dat$Abundance[rel_abund_dat$species=='target species'&
                          rel_abund_dat$mean_attract=='linear'] <-
  (c(300, 400, 500, 600, 700, 800))[rel_abund_dat$Year[rel_abund_dat$species=='target species'&
                       rel_abund_dat$mean_attract=='linear']]

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
                 constant = 'constant target species \nabundance',
                 `linear` = 'linearly increasing slow \nspecies abundance'
               ),
               sat_level = c(
                 high = 'saturation "common"'
               ))) +
  ylab('Mean absolute abundance')  +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  guides(fill='none') +
  theme(strip.background = element_blank(),strip.text = element_blank(),
        legend.position = c(0.43,0.94),legend.box.background=element_blank(),
        legend.background=element_blank()) +
  guides(linetype=guide_legend('')) +
  ylim(c(0,800))

# THESE ARE DESIGNED FOR A4 LANDSCAPE or 5.83 x 11.3
# NOTICE THE HACK IN MULTIPLOT'S LAYOUT ARGUMENT

# First show the 'typical' performance of methods
multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
Results %>%
  group_by(sat_effects ,mean_attract, nsim) %>%
  mutate(all_converged = mean(Converge)) %>%
  ungroup() %>%
  filter(Station>1, all_converged==1) %>%
  group_by(model, Station, sat_effects, mean_attract) %>%
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
  facet_grid(mean_attract ~sat_effects , scales = 'free_y',
             labeller = labeller(
               sat_effects=c(
                 `no saturation` = 'p* = 1',
                 saturation = 'p* = 0.85'
               ),
               bite_fun=c(
                 constant = 'Uniform Distributions',
                 mixed = 'Mixture of Distributions'
               ),
               mean_attract = c(
                 constant = 'constant target species \nabundance',
                 linear = 'linearly increasing target \nspecies abundance '
               ),
               sat_level = c(
                 low = 'saturation less "common"',
                 high = 'saturation "common"'
               )
             )) +
  geom_hline(yintercept=0) +
  ggtitle('Bias in relative abundance for each method')+#,
          #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe hook location abilities after 85% of baits removed.') +
  ylab('Bias') +
  xlab('Year') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
        legend.position = 'left', strip.text.y = element_blank()) +
  scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q')) +
  scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q'),
                     values=c('circle','triangle','square','square')) +
  guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
         shape=guide_legend(title = 'Method')),
layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    group_by(sat_effects ,mean_attract, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, sat_effects, mean_attract) %>%
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
    facet_grid(mean_attract ~sat_effects , scales = 'free_y',
               labeller = labeller(
                 sat_effects=c(
                   `no saturation` = 'p* = 1',
                   saturation = 'p* = 0.85'
                 ),
                 bite_fun=c(
                   constant = 'Uniform Distributions',
                   mixed = 'Mixture of Distributions'
                 ),
                 mean_attract = c(
                   constant = 'constant target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0) +
    ggtitle('Bias in relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe hook location abilities after 85% of baits removed.') +
    ylab('Bias') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title = 'Method')),
  layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    group_by(sat_effects ,mean_attract, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, sat_effects, mean_attract) %>%
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
    facet_grid(mean_attract ~sat_effects , scales = 'free_y',
               labeller = labeller(
                 sat_effects=c(
                   `no saturation` = 'p* = 1',
                   saturation = 'p* = 0.85'
                 ),
                 bite_fun=c(
                   constant = 'Uniform Distributions',
                   mixed = 'Mixture of Distributions'
                 ),
                 mean_attract = c(
                   constant = 'constant target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0.95) +
    ylim(c(0,1.05)) +
    ggtitle('Coverage of intervals of relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe hook location abilities after 85% of baits removed.') +
    ylab('Coverage') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title = 'Method')),
  layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))


multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    group_by(sat_effects ,mean_attract, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, sat_effects, mean_attract) %>%
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
    facet_grid(mean_attract ~sat_effects , scales = 'free_y',
               labeller = labeller(
                 sat_effects=c(
                   `no saturation` = 'p* = 1',
                   saturation = 'p* = 0.85'
                 ),
                 bite_fun=c(
                   constant = 'Uniform Distributions',
                   mixed = 'Mixture of Distributions'
                 ),
                 mean_attract = c(
                   constant = 'constant target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0) +
    ggtitle('MSE in relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe hook location abilities after 85% of baits removed.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title = 'Method')),
  layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    group_by(sat_level, mean_attract, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1, model!='naive') %>%
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
                   constant = 'Uniform Distributions',
                   mixed = 'Mixture of Distributions'
                 ),
                 mean_attract = c(
                   constant = 'constant target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    geom_hline(yintercept=0) +
    ggtitle('MSE in relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe hook location abilities after 85% of baits removed.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q'),
                       values=c('circle','triangle','square','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title = 'Method')),
  layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))


multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    group_by(sat_level, mean_attract, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
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
                   constant = 'Uniform Distributions',
                   mixed = 'Mixture of Distributions'
                 ),
                 mean_attract = c(
                   constant = 'constant target species \nabundance',
                   linear = 'linearly increasing target \nspecies abundance '
                 ),
                 sat_level = c(
                   low = 'saturation less "common"',
                   high = 'saturation "common"'
                 )
               )) +
    ggtitle('Proportion of Fishing Events with Specified Level of Hook Saturation')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation.\nColumns describe arrival time distributions and hook location abilities after 85% of baits removed.') +
    ylab('Proportion of fishing events') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('100%','>85%')) +
    guides(color=guide_legend(override.aes=list(fill=NA), title = 'Method'),
           shape=guide_legend(title = 'Method')),
    labs(colour='Baits Removed'),
  layout = matrix(c(rep(NA,40),rep(1,460),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

