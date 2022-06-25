# Simulation study demonstrating the upper bound used
# in the sim experiments is conservative
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)
seed <- 25042021 # date
set.seed(seed)
nspecies <- 3
# bite_funs <- cbind(rep('constant',3),
#                    c('exp_decay','mixed', 'constant'))
bite_funs <-c('exp_decay','mixed', 'constant')
soak_time <- 5
n_hooks <- 800
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
n_sim <- 1
hook_sat_level <- 0.85 # true proportion at which saturation effects begin
cprop=0.85 # assumed proportion at which saturation effects begin
#cprop2=0.95 # assumed proportion "" "" second model
upper_bound_quantiles <- c(0.85, 0.95, 1)
#lower_bound_quantile <- c(0.1)

# how much does the bite rate of each species (linearly) decrease at 100% saturation
# assume linear decrease from 85% saturation onwards
# saturation_effects <- cbind(c(0, 0, 0),
#                            c(0, 0.2, 0.8))
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

# plot(x=seq(from=0,to=1, length.out=100), y=sat_fun(0.8,seq(from=0,to=1, length.out=100)))
#
# # Plot the bite rates of the two species
# plot(x=seq(from=0, to=5, length.out=100),y=exp(-seq(from=0, to=5, length.out=100))/(1-exp(-5)),
#      type = 'l', main='bite rate functions', ylab='bite rate',xlab='time')
# lines(x=seq(from=0, to=5, length.out=100),y=rep(0.2, times=100), col='red')

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

# define function for competition
# competition_fun <-
#   function(data, p_sat, p_baseline){
#     val <- rep(1, length(p_sat))
#
#     preddata=data
#     preddata$prop_sat=p_sat
#     #preddata$region <- preddata$region_INLA
#     preddata2=preddata
#     preddata2$prop_sat=p_baseline
#     #browser()
#     val2 <-
#       apply(predict.gam(mod_species1,
#                         newdata=preddata,
#                         type = 'terms'),
#             1, FUN = function(x){exp(sum(x))})/
#       apply(predict.gam(mod_species1,
#                         newdata=preddata2,
#                         type = 'terms'),
#             1, FUN = function(x){exp(sum(x))})
#     # predict.gam(mod_species1,
#     #             newdata=preddata,
#     #             type = 'response') /
#     # predict.gam(mod_species1,
#     #             newdata=preddata2,
#     #             type = 'response')
#
#     val[which(p_sat > p_baseline)] <-
#       val2[which(p_sat > p_baseline)]
#     return(val)
#   }

Results <- data.frame(
  nsim = rep(1:n_sim, each=nstation*2*2*4*nyears),
  sat_level = rep(rep(c('low','high'), times=n_sim),each=2*4*nstation*nyears),
  mean_attract = rep(rep(mean_attract, times=n_sim*2), each=4*nstation*nyears),
  correlation = rep(rep(c('negative','low','medium','high'),times=n_sim*2*2), each=nstation*nyears),
  #sat_effects = rep(rep(c('no saturation','saturation'),times=n_sim*2*2), each=2*4*nstation),
  #bite_fun = rep(rep(c('constant','mixed'), times=n_sim*2*2*2), each=4*nstation),
  #model=rep(rep(c('naive','adjust','censored_upper85','censored_upper95','censored_upper100','censored'), times=n_sim*2*2*4), each=nstation*nyears),
  N_Attracted=rep(0, times=nyears*n_sim*2*2*4*nstation),
  Upper_Bound=rep(0, times=nyears*n_sim*2*2*4*nstation),
  Censored=rep(0, times=nyears*n_sim*2*2*4*nstation),
  Prop_Sat=rep(0, times=nyears*n_sim*2*2*4*nstation)
  )

results_mapper <- function(n,i,j,k,mod)
{
  return(
    which(Results$nsim == n & Results$sat_level == c('low','high')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$correlation == c('negative','low','medium','high')[k])
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
        # NOTE THAT WE NOW SCALE THE MEAN OF THE TARGET SPECIES BY THE EXPONENTIAL OF THE
        # LOG-NORMAL VARIANCE TO ENSURE THE SAME MEAN FROM THE PREVIOUS EXPERIMENT
        # DOESN'T AFFECT relative abundance!!
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
          # for(l in 1:dim(bite_funs)[2])
          # {
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
                                     #exp(rnorm(dim(nbite)[1], mean = 0, sd=sd_log_bite[nbite$species])))

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

          # Use the Baranov Catch equation to (starting at 85% saturation) to derive scale factor
          dat <- nbite[nbite$species==3,]
          scale_fac <- rep(0, length(dat$bites))
          scale_fac[dat$prop_sat>cprop] <-
            comp_factor_fun(1-signif((dat[dat$prop_sat>cprop,]$prop_sat-cprop)/(1-cprop),5),
                            rep(round((1-cprop)*n_hooks),sum(dat$prop_sat>cprop)))

          upper_bound <- rep(0,length(dat$bites))
          upper_bound[dat$prop_sat>cprop] <- round(
            (dat$prop_sat[dat$prop_sat>cprop]-cprop)*n_hooks*
              scale_fac[dat$prop_sat>cprop])

          Results[results_mapper(nsim,i,j,k),'N_Attracted'] <- as.numeric(dat$attracted)
          Results[results_mapper(nsim,i,j,k),'Upper_Bound'] <- dat$bites+upper_bound
          Results[results_mapper(nsim,i,j,k),'Censored'] <- ifelse(dat$prop_sat>cprop,1,0)
          Results[results_mapper(nsim,i,j,k),'Prop_Sat'] <- dat$prop_sat

      }
    }
  }
}

# Change the factors to ordered factors to improve plotting
Results$sat_level <- factor(Results$sat_level, levels=c('low','high'), ordered = T)
Results$mean_attract <- factor(Results$mean_attract, levels=c('constant','linear'), ordered = T)
Results$correlation <- factor(Results$correlation, levels=c('negative','low','medium','high'), ordered = T)

Results %>%
  group_by(sat_level, mean_attract, correlation, nsim) %>%
  filter(Censored==1) %>%
  ggplot(aes(x=N_Attracted, y=Upper_Bound)) +
  geom_point(position = position_jitter(width=1, height = 1)) +
  geom_abline(slope=1, intercept = 0) +
  facet_grid(correlation~sat_level+mean_attract)
# Success! The vast majority of the upper bound estimates lie
# far above the true (but unobserved) number of attracted individuals

Results %>%
  group_by(sat_level, mean_attract, correlation, nsim) %>%
  filter(Censored==1) %>%
  ggplot(aes(x=Prop_Sat, y=Upper_Bound-N_Attracted)) +
  geom_point(position = position_jitter(width=0, height = 1)) +
  geom_abline(slope=0, intercept = 0) +
  facet_grid(correlation~sat_level+mean_attract)

# Precisely what fraction lie above and by what amount?
Results %>%
  group_by(sat_level, mean_attract, correlation, nsim) %>%
  filter(Censored==1) %>%
  summarize(Mean_Amount_Above = mean(Upper_Bound-N_Attracted),
            Proportion_Above = mean((Upper_Bound-N_Attracted)>0)) %>%
  View()
# Great - across all settings, the upper bound > n_attracted in vast majority
# of times. The mean amount above is >250

Results %>%
  filter(Censored==1) %>%
  summarize(Mean_Amount_Above = mean(Upper_Bound-N_Attracted),
            Proportion_Above = mean((Upper_Bound-N_Attracted)>0))
# 99.9% of time - mean amount 350
