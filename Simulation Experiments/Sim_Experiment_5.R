# Simulation study demonstrating the censored hook competition method
# 3 species, as in experiment 3, but the target species reaches hook first
# We comment this code less due to the overlap with sim experiment 2's code
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)
seed <- 25032022 # date
set.seed(seed)
nspecies <- 3
# Note that the target species (3rd) is now fastest to arrive to baits with exp dist
bite_funs <- c('constant', 'mixed', 'exp_decay')
soak_time <- 5
n_hooks <- 800
# definition of stations and years are opposite to the manuscript
nstation <- 6
nyears <- 100
saturation_level <- c(1,2)
mean_attract <- c('constant', 'linear')
cov_matrices <- list( neg=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,-0.6,0,-0.6,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2)),
                      low=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,0.1,0,0.1,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2)),
                      med=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,0.3,0,0.3,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2)),
                      high=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,0.6,0,0.6,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2)))
n_sim <- 100
hook_sat_level <- 0.85 # true proportion at which saturation effects begin
cprop=0.95 # assumed proportion at which saturation effects begin
upper_bound_quantiles <- c(0.85, 0.95, 1)

# create saturation function
sat_fun <- function(sat_effect, sat_level, hook_sat_level=0.85)
{
  val <- rep(1, length(sat_level))

  val[which(sat_level>hook_sat_level)] <-
    (1 - sat_effect*((sat_level[which(sat_level>hook_sat_level)]-hook_sat_level)/(1-hook_sat_level)))
  return(val)
}

# try the hook competition ICR method
comp_factor_fun <- function(prop_hook, n_hook)
{
  prop <- prop_hook
  # if all hooks saturated - map to 1 hook
  prop[which(prop == 0)] <- 1 / n_hook[which(prop == 0)]
  return(-log(prop)/(1-prop))
}

# sample the unICR bite times for each of the 'attracted' fish for each fishing event
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
  nsim = rep(1:n_sim, each=nstation*2*2*4*6),
  sat_level = rep(rep(c('low','high'), times=n_sim),each=2*4*6*nstation),
  mean_attract = rep(rep(mean_attract, times=n_sim*2), each=4*6*nstation),
  correlation = rep(rep(c('negative','low','medium','high'),times=n_sim*2*2), each=6*nstation),
  model=rep(rep(c('CPUE','adjust','censored_upper85','censored','censored_upper85_cprop1','censored_cprop1'), times=n_sim*2*2*4), each=nstation),
  Bias=rep(0, times=n_sim*2*2*4*6*nstation),
  Converge=rep(0, times=n_sim*2*2*4*6*nstation),
  RMSE=rep(0, times=n_sim*2*2*4*6*nstation),
  Coverage=rep(0, times=n_sim*2*2*4*6*nstation),
  Rel_Bias=rep(0, times=n_sim*2*2*4*6*nstation),
  Rel_RMSE=rep(0, times=n_sim*2*2*4*6*nstation),
  Rel_Coverage=rep(0, times=n_sim*2*2*4*6*nstation),
  Station=rep(1:nstation, times=n_sim*2*2*4*6),
  Prop_Sat_85=rep(0, times=n_sim*2*2*4*6*nstation),
  Prop_Sat_100=rep(0, times=n_sim*2*2*4*6*nstation))

results_mapper <- function(n,i,j,k,mod)
{
  return(
    which(Results$nsim == n & Results$sat_level == c('low','high')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$correlation == c('negative','low','medium','high')[k] &
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
      for(k in 1:length(cov_matrices))
      {
        if(mean_attract[j] == 'constant')
        {
          mean_bite_gen =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (rep(100, 6)/exp(cov_matrices[[k]][3,3]/2)))
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (rep(100, 6)))
        }
        if(mean_attract[j] == 'linear')
        {
          mean_bite_gen =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (c(120, 140, 160, 180, 200, 220)-100)/exp(cov_matrices[[k]][3,3]/2))
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (c(120, 140, 160, 180, 200, 220)-100))
        }
        saturation_effect <- c(0, 0.2, 0.8)
        # sample the number of each species that WOULD bite at each station for each year if hooks were available
        nbite <- data.frame(bites = rep(0, times=nspecies*nstation*nyears),
                            attracted = rep(0, times=nspecies*nstation*nyears),
                            species=rep(1:3, each=nstation*nyears),
                            station=rep(1:nstation, times=nspecies*nyears),
                            year=rep(rep(1:nyears, each=nstation), times=nspecies))
        nbite$attracted <- rpois(dim(nbite)[1],
                                 lambda = as.numeric(mean_bite_gen[cbind(nbite$station,nbite$species)])*
                                   exp(as.numeric(rmvn(n=nstation*nyears, mu = rep(0,nspecies), V=cov_matrices[[k]])[cbind(rep(1:(nstation*nyears),times=nspecies),nbite$species)])))
        for(i2 in 1:nstation)
        {
          for(j2 in 1:nyears)
          {
            bite_time_1 <- bite_samp(bite_funs[1],sum(nbite$attracted[nbite$species==1 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_1)>soak_time)
            {
              bite_time_1[bite_time_1>soak_time] <-
                bite_samp(bite_funs[1],sum(bite_time_1>soak_time))
            }
            # repeat for species 2 from uniform distribution
            bite_time_2 <- bite_samp(bite_funs[2],sum(nbite$attracted[nbite$species==2 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_2)>soak_time)
            {
              bite_time_2[bite_time_2>soak_time] <-
                bite_samp(bite_funs[2],sum(bite_time_2>soak_time))
            }

            # repeat for species 3 from uniform distribution
            bite_time_3 <- bite_samp(bite_funs[3],sum(nbite$attracted[nbite$species==3 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_3)>soak_time)
            {
              bite_time_3[bite_time_3>soak_time] <-
                bite_samp(bite_funs[3],sum(bite_time_3>soak_time))
            }

            # Now we sample the first n_hooks*n_hook_sat_level unICR
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
              # Now we sample the first n_hooks*n_hook_sat_level unICR
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

        Results[results_mapper(nsim,i,j,k,'CPUE'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
        Results[results_mapper(nsim,i,j,k,'adjust'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
        Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
        Results[results_mapper(nsim,i,j,k,'CPUE'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'adjust'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
        Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))
        Results[results_mapper(nsim,i,j,k,'censored'),'Prop_Sat_100'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat==1)}))
        Results[results_mapper(nsim,i,j,k,'censored'),'Prop_Sat_85'] <- as.numeric(by(nbite, nbite$station, FUN=function(x){mean(x$prop_sat>cprop)}))

        # fit a CPUE model that ignores competition
        dat <- nbite[nbite$species == 3,]
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
        if(sum(is.na(mod2$summary.fixed$mean))>0 | sum(mod2$misc$mode.status!=0)>0 )
        {
          mod2 <- NULL
        }

        if(!is.null(mod2))
        {
          Results[results_mapper(nsim,i,j,k,'CPUE'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod2,
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

          Results[results_mapper(nsim,i,j,k,'CPUE'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'CPUE'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'CPUE'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'CPUE'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'CPUE'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'CPUE'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        # ICR method
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
          Results[results_mapper(nsim,i,j,k,'adjust'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod4,
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

          Results[results_mapper(nsim,i,j,k,'adjust'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'adjust'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'adjust'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'adjust'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'adjust'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'adjust'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        # censored method
        upper_bound <- rep(0, length(nbite[nbite$species==3,]$prop_sat))
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

        }
        dat <- nbite[nbite$species==3,]
        dat$low <- rep(Inf,dim(dat)[1])
        dat$high <- rep(Inf,dim(dat)[1])

        quant_regions <-
          as.numeric(by(dat$bites,
                        dat$station,
                        FUN = function(x){quantile(x,0.85, na.rm=T)}, simplify = T))

        quant_regions2 <- rep(0, nstation)

        dat$low[which(dat$prop_sat >= cprop &
                        0 < upper_bound &
                        dat$bites <= quant_regions[dat$station] &
                        dat$bites >= quant_regions2[dat$station])] <-
          as.matrix(dat[which(dat$prop_sat >= cprop &
                                0 < upper_bound &
                                dat$bites <= quant_regions[dat$station] &
                                dat$bites >= quant_regions2[dat$station]),
                        c('bites')])[,1]

        dat$high[which(dat$prop_sat >= cprop &
                         0 < upper_bound &
                         dat$bites <= quant_regions[dat$station] &
                         dat$bites >= quant_regions2[dat$station])] <-
          dat$bites[which(dat$prop_sat >= cprop &
                            0 < upper_bound &
                            dat$bites <= quant_regions[dat$station] &
                            dat$bites >= quant_regions2[dat$station])] +
          upper_bound[which(dat$prop_sat >= cprop &
                              0 < upper_bound &
                              dat$bites <= quant_regions[dat$station] &
                              dat$bites >= quant_regions2[dat$station])]


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
          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod6,
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

          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'censored_upper85'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        # Repeat for other censored models with differing upper observations
        dat <- nbite[nbite$species==3,]
        dat$low <- rep(Inf,dim(dat)[1])
        dat$high <- rep(Inf,dim(dat)[1])

        quant_regions2 <- rep(0, nstation)

        dat$low[which(dat$prop_sat >= 1 &
                        dat$bites <= quant_regions[dat$station] &
                        dat$bites >= quant_regions2[dat$station])] <-
          as.matrix(dat[which(dat$prop_sat >= 1 &
                                dat$bites <= quant_regions[dat$station] &
                                dat$bites >= quant_regions2[dat$station]),
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
          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod8,
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

          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        # Another censored method
        dat <- nbite[nbite$species==3,]
        dat$low <- rep(Inf,dim(dat)[1])
        dat$high <- rep(Inf,dim(dat)[1])

        quant_regions2 <- rep(0, nstation)

        dat$low[which(dat$prop_sat == 1)] <-
          as.matrix(dat[which(dat$prop_sat ==1),
                        c('bites')])[,1]

        ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

        dat$event_ID <- 1:dim(dat)[1]

        mod10 <-
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
        if(sum(is.na(mod10$summary.fixed$mean))>0 | sum(mod10$misc$mode.status!=0)>0 )
        {
          mod10 <- NULL
        }

        if(!is.null(mod10))
        {
          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod10,
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

          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'censored_cprop1'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        # Another censored method
        dat <- nbite[nbite$species==3,]
        dat$low <- rep(Inf,dim(dat)[1])
        dat$high <- rep(Inf,dim(dat)[1])

        dat$low[which(dat$prop_sat >= cprop &
                        0 < upper_bound)] <-
          as.matrix(dat[which(dat$prop_sat >= cprop &
                                0 < upper_bound),
                        c('bites')])[,1]

        dat$high[which(dat$prop_sat >= cprop &
                         0 < upper_bound)] <-
          dat$bites[which(dat$prop_sat >= cprop &
                            0 < upper_bound)] +
          upper_bound[which(dat$prop_sat >= cprop &
                              0 < upper_bound)]


        ind_resp <- which(names(dat) %in% c('bites', 'low', 'high'))

        dat$event_ID <- 1:dim(dat)[1]


        mod12 <-
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
        if(sum(is.na(mod12$summary.fixed$mean))>0 | sum(mod12$misc$mode.status!=0)>0 )
        {
          mod12 <- NULL
        }

        if(!is.null(mod12))
        {
          Results[results_mapper(nsim,i,j,k,'censored'),'Converge'] <- 1

          parameters <- inla.posterior.sample(2000,mod12,
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

          Results[results_mapper(nsim,i,j,k,'censored'),'Rel_Bias'][-1] <-
            parameters[2:6,2] -
            mean_bite[-1,3]/mean_bite[1,3]

          Results[results_mapper(nsim,i,j,k,'censored'),'Rel_RMSE'][-1] <-
            (parameters[2:6,2] -
               mean_bite[-1,3]/mean_bite[1,3])^2

          Results[results_mapper(nsim,i,j,k,'censored'),'Rel_Coverage'][-1] <-
            ifelse(parameters[2:6,1] <= mean_bite[-1,3]/mean_bite[1,3] & parameters[2:6,3] >= mean_bite[-1,3]/mean_bite[1,3],
                   1,0)
          Results[results_mapper(nsim,i,j,k,'censored'),'Bias'] <-
            parameters[c(1,7:11),2] - mean_bite_gen[,3]
          Results[results_mapper(nsim,i,j,k,'censored'),'RMSE'] <-
            (parameters[c(1,7:11),2] - mean_bite_gen[,3])^2
          Results[results_mapper(nsim,i,j,k,'censored'),'Coverage'] <-
            ifelse(parameters[c(1,7:11),1] <= mean_bite[,3] & parameters[c(1,7:11),3] >= mean_bite[,3],
                   1,0)
        }

        print(Results[results_mapper(nsim,i,j,k,'CPUE'),])
        print(Results[results_mapper(nsim,i,j,k,'adjust'),])
        print(Results[results_mapper(nsim,i,j,k,'censored_upper85'),])
        print(Results[results_mapper(nsim,i,j,k,'censored_upper85_cprop1'),])
        print(Results[results_mapper(nsim,i,j,k,'censored_cprop1'),])
        print(Results[results_mapper(nsim,i,j,k,'censored'),])
        rm(mod6,mod8,mod10,mod12)
      }
    }
  }
}
# Save results
#saveRDS(Results, 'Simulation_Results_Correlated_Corrected_Aggressive.rds')
Results <- readRDS('Simulation_Results_Correlated_Corrected_Aggressive.rds')

library(inlabru)
# Create artificial 'relative abundance' of target and slow species plots
rel_abund_dat <- data.frame(expand.grid(
  species=c('target species','slow species'),
  sat_level=factor(c('low','high'), levels=c('low','high'), ordered = T),
  mean_attract=factor(c('constant','linear')),
  Year=c(1,2,3,4,5,6)))
rel_abund_dat$Abundance <- 1
rel_abund_dat$Abundance[rel_abund_dat$species=='target species'&
                          rel_abund_dat$mean_attract=='linear'] <-
  rel_abund_dat$Year[rel_abund_dat$species=='target species'&
                       rel_abund_dat$mean_attract=='linear']
rel_abund_dat$Abundance[rel_abund_dat$species=='slow species'&
                          rel_abund_dat$sat_level=='low'] <-
  (c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='slow species'&
                                                            rel_abund_dat$sat_level=='low']]
rel_abund_dat$Abundance[rel_abund_dat$species=='slow species'&
                          rel_abund_dat$sat_level=='high'] <-
  2*(c(120,140, 160, 180, 200, 280)/120)[rel_abund_dat$Year[rel_abund_dat$species=='slow species'&
                                                              rel_abund_dat$sat_level=='high']]

rel_abund_plot <-
  ggplot(rel_abund_dat,
         aes(x=Year, y=Abundance, linetype=species) ) +
  geom_rect(data= ~.x[.x$Year==1,],
            aes(x=Year, y=Abundance, linetype=species, fill = mean_attract),
            xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_line() +
  scale_fill_discrete(labels=c('Target Species', 'Slow Species')) +
  facet_grid(mean_attract+sat_level~.,
             labeller = labeller(
               mean_attract = c(
                 constant = 'constant target species \nabundance',
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

# Change the factors to ordered factors to improve plotting
Results$sat_level <- factor(Results$sat_level, levels=c('low','high'), ordered = T)
Results$mean_attract <- factor(Results$mean_attract, levels=c('constant','linear'), ordered = T)
Results$correlation <- factor(Results$correlation, levels=c('negative','low','medium','high'), ordered = T)
Results$model <- factor(Results$model, levels=c('CPUE','adjust','censored','censored_cprop1','censored_upper85_cprop1','censored_upper85'), ordered = T)

# THESE ARE DESIGNED FOR 5.83x11.3
# NOTICE THE HACK IN MULTIPLOT'S LAYOUT ARGUMENT
multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
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
    #geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation. Columns are correlation \nof target species\' abundance with saturation events. 95% Confidence intervals smaller than plotting shapes so omitted.') +
    ylab('Bias') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','triangle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
    mutate(Mean = median(Rel_Bias),
           UCL = quantile(Rel_Bias, probs=0.975),
           LCL = quantile(Rel_Bias, probs=0.025)) %>%
    ungroup(Station) %>%
    mutate(random_samp = c(T, rep(F, length(Rel_Bias)-1))) %>%
    group_by(Station) %>%
    ggplot(aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, shape=model)) +
    geom_rect(data= ~.x[.x$random_samp==T,],
              aes(x=Station, y=Mean, ymin=LCL, ymax=UCL, colour=model, group=model, fill = mean_attract),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3) +
    #geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation. Columns are correlation \nof target species\' abundance with saturation events. 95% Confidence intervals smaller than plotting shapes so omitted.') +
    ylab('Bias') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','triangle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored','adjust'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
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
    #geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation. Columns are correlation \nof target species\' abundance with saturation events. 95% Confidence intervals smaller than plotting shapes so omitted.') +
    ylab('Bias') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
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
    geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
    ggtitle('Coverage of intervals of relative abundance for each method')+#,
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation,\nColumns are correlation of target species\' abundance with saturation events') +
    ylab('Coverage') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','triangle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored','adjust'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
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
    #geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation. Columns are correlation \nof target species\' abundance with saturation events. 95% Confidence intervals smaller than plotting shapes so omitted.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))

multiplot(
  rel_abund_plot + ggtitle('Simulated abundance'),
  Results %>%
    filter(!(model %in% c('censored'))) %>%
    group_by(sat_level, mean_attract, correlation, nsim) %>%
    mutate(all_converged = mean(Converge)) %>%
    ungroup() %>%
    filter(Station>1, all_converged==1) %>%
    group_by(model, Station, correlation, sat_level,  mean_attract) %>%
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
    #geom_errorbar(position = position_dodge(width=0.9)) +
    geom_point(position = position_dodge(width=0.9), size=2) +
    facet_grid(mean_attract + sat_level ~ correlation , scales = 'free_y',
               labeller = labeller(
                 correlation=c(
                   negative = 'negative correlation',
                   low = 'low positive correlation',
                   medium = 'med positive correlation',
                   high = 'high positive correlation'
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
            #subtitle = 'Rows are trends in relative abundance of target species and the average degree of hook saturation. Columns are correlation \nof target species\' abundance with saturation events. 95% Confidence intervals smaller than plotting shapes so omitted.') +
    ylab('MSE') +
    xlab('Year') + guides(fill='none') +
    scale_fill_brewer(palette = 'Pastel1') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          legend.position = 'left', strip.text.y = element_blank()) +
    scale_color_viridis_d(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q')) +
    scale_shape_manual(labels=c('CPUE','ICR','Censored 1','Censored 1 Q','Censored 0.95 Q'),
                       values=c('circle','triangle','cross','plus','square')) +
    guides(color=guide_legend(override.aes=list(fill=NA)))+ labs(colour='Method', shape='Method'),
  layout = matrix(c(rep(NA,32),rep(1,468),rep(2,2000)), nrow = 500, ncol = 5, byrow = F))
