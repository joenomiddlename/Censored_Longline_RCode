# Simulation study to identify a 'signature' for choosing hatp*
# We simulate catch counts using data generation steps similar to
# experiments 3 and 5
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)
seed <- 17042022#25032022 # date
set.seed(seed)
nspecies <- 3
# This time we allow the target species to be fast, equally matched, and slow to reach baits
bite_funs <- list(aggressive=c('constant','mixed', 'exp_decay'),
                  equally_matched=rep('constant',3),
                  slow=c('exp_decay','mixed', 'constant')
)

soak_time <- 5
n_hooks <- 800
# Note that the definitions of station and year are opposite to manuscript
nstation <- 6
nyears <- 100
# Below: either the 'other' species group or the target species is highly abundant
saturation_level <- cbind(c(1,2,1),c(1,1,2))
# What are the temporal abundance trends of the 3 species.
# 'constant' means all three species have constant abundance, etc.,
mean_attract <- c('constant', 'linear target', 'linear aggressive')
# Species are uncorrelated but overdispersed.
sd_species <- list(equal = rep(0.2,3),
                   overdisp_target=c(0.2, 0.2, 0.8),
                   overdisp_other=c(0.8, 0.2, 0.2))
# Only 20 simulation iterations done here
n_sim <- 20
hook_sat_level <- 0.85# true proportion at which saturation effects begin
cprop <- 0.95 # The value to use to test for breakdown
# how strong are the saturation effects for each species? If 0, then no hook competition effects
saturation_effect <- list(target=c(0, 0.2, 0.8),
                          none=rep(0,3))

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
  nsim = rep(1:n_sim, each=3*3*3*3*3*2),
  sat_level = rep(rep(c('low','high other','high target'), times=n_sim),each=3*3*3*3*2),
  mean_attract = rep(rep(mean_attract, times=n_sim*3), each=3*3*3*2),
  correlation = rep(rep(c('negative','none','positive'),times=n_sim*3*3), each=3*3*2),
  target_species=rep(rep(c('aggressive','equally_matched','slow'),times=n_sim*3*3*3), each=3*2),
  sd_species=rep(rep(c('equal','overdisp_target','overdisp_other'),times=n_sim*3*3*3*3), each=2),
  sat_effects = rep(c('target','no saturation'),times=n_sim*3*3*3*3),
  Delta_Mean_Catch_Count=0
)

results_mapper <- function(sim,i,j,k,l,m,n)
{
  return(
    which(  Results$nsim == sim &
            Results$sat_level == c('low','high other','high target')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$correlation == c('negative','none','positive')[k] &
            Results$target_species == c('aggressive','equally_matched','slow')[l] &
            Results$sd_species == c('equal','overdisp_target','overdisp_other')[m] &
            Results$sat_effects == c('target','no saturation')[n])
  )
}

for(sim in 1:n_sim)
{
  print(paste0('simulation number ',sim, ' out of ',n_sim, ' processing...'))
for(i in 1:dim(saturation_level)[1])
{
  for(j in 1:length(mean_attract))
  {
    for(k in 1:length(unique(Results$correlation)))
    {
      for(l in 1:length(unique(Results$target_species)))
      {
        for(m in 1:length(unique(Results$sd_species)))
        {
          for(n in 1:length(unique(Results$sat_effects)))
          {


            cov_matrices <- list( neg=
                                    diag(sd_species[[m]]) %*%
                                    matrix(c(1,0,0,0,1,-0.6,0,-0.6,1), 3,3,byrow = T) %*%
                                    diag(sd_species[[m]]),
                                  none=
                                    diag(sd_species[[m]]) %*%
                                    matrix(c(1,0,0,0,1,0,0,0,1), 3,3,byrow = T) %*%
                                    diag(sd_species[[m]]),
                                  positive=
                                    diag(sd_species[[m]]) %*%
                                    matrix(c(1,0,0,0,1,0.6,0,0.6,1), 3,3,byrow = T) %*%
                                    diag(sd_species[[m]]))

            if(mean_attract[j] == 'constant')
            {
              mean_bite_gen =
                cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i,1],
                      rep(400,6),
                      (rep(100, 6)/exp(cov_matrices[[k]][3,3]/2))*saturation_level[i,2])
              mean_bite =
                cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i,1],
                      rep(400,6),
                      (rep(100, 6))*saturation_level[i,2])
            }
            if(mean_attract[j] == 'linear target')
            {
              mean_bite_gen =
                cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i,1],
                      rep(400,6),
                      ((c(120, 140, 160, 180, 200, 220)-100)/exp(cov_matrices[[k]][3,3]/2))*saturation_level[i,2])
              mean_bite =
                cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i,1],
                      rep(400,6),
                      (c(120, 140, 160, 180, 200, 220)-100)*saturation_level[i,2])
            }
            if(mean_attract[j] == 'linear aggressive')
            {
              mean_bite_gen =
                cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i,1],
                      rep(400,6),
                      rep(100, 6)*saturation_level[i,2])
              mean_bite =
                cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i,1],
                      rep(400,6),
                      rep(100, 6)*saturation_level[i,2])
            }

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
                bite_time_1 <- bite_samp(bite_funs[[l]][1],sum(nbite$attracted[nbite$species==1 &
                                                                                 nbite$station==i2 &
                                                                                 nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_1)>soak_time)
                {
                  bite_time_1[bite_time_1>soak_time] <-
                    bite_samp(bite_funs[[l]][1],sum(bite_time_1>soak_time))
                }
                # repeat for species 2 from uniform distribution
                bite_time_2 <- bite_samp(bite_funs[[l]][2],sum(nbite$attracted[nbite$species==2 &
                                                                                 nbite$station==i2 &
                                                                                 nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_2)>soak_time)
                {
                  bite_time_2[bite_time_2>soak_time] <-
                    bite_samp(bite_funs[[l]][2],sum(bite_time_2>soak_time))
                }

                # repeat for species 3 from uniform distribution
                bite_time_3 <- bite_samp(bite_funs[[l]][3],sum(nbite$attracted[nbite$species==3 &
                                                                                 nbite$station==i2 &
                                                                                 nbite$year==j2]))
                # truncate them to 0-5 interval
                while(max(bite_time_3)>soak_time)
                {
                  bite_time_3[bite_time_3>soak_time] <-
                    bite_samp(bite_funs[[l]][3],sum(bite_time_3>soak_time))
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
                                                    prob=sat_fun(sat_effect = saturation_effect[[n]][1],
                                                                 sat_level = current_sat_level,
                                                                 hook_sat_level = hook_sat_level))==1,
                                             time_ind[k2], NA)
                      flag <- F
                    }
                    if(species_ind[time_ind[k2]]==2 & flag)
                    {
                      time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                    prob=sat_fun(sat_effect = saturation_effect[[n]][2],
                                                                 sat_level = current_sat_level,
                                                                 hook_sat_level = hook_sat_level))==1,
                                             time_ind[k2], NA)
                      flag <- F
                    }
                    if(species_ind[time_ind[k2]]==3 & flag)
                    {
                      time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                    prob=sat_fun(sat_effect = saturation_effect[[n]][3],
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

            # NEW CODE

            # The model is a regression to see if the mean
            # catch counts decreases/increases at high values of gear saturation

            # Note that High_Sat is defined as p_it > 0.95 and Medium Sat as 0.8 < p_it <= 0.95
            mod2 <-
              inla(data=
                     nbite %>%
                     filter(species==3) %>%
                     mutate(event_ID = 1:row_number(),
                            High_Sat=ifelse(prop_sat > cprop,1,0),
                            Medium_Sat=ifelse(prop_sat > cprop-0.15 & prop_sat <= cprop,1,0)),
                   formula = bites ~ High_Sat + Medium_Sat + factor(station) + f(event_ID,model='iid'),
                   family='poisson',
                   control.compute = list(config=TRUE),
                   control.mode = list(theta=5, restart=T)
              )
            summary(mod2)
            # What is the mean change in catch count from medium to high saturation
            test <-
              apply(
                inla.posterior.sample.eval(fun=function(x){High_Sat-Medium_Sat},
                                           samples=inla.posterior.sample(n=500, mod2, selection = list('High_Sat'=1,'Medium_Sat'=1))),
                1, FUN = function(x){quantile(x,probs=0.5)}
                )

            Results[results_mapper(sim,i,j,k,l,m,n),'prop_sat'] <- mean(nbite$prop_sat)
            Results[results_mapper(sim,i,j,k,l,m,n),'Delta_Mean_Catch_Count'] <- test

            print(Results[results_mapper(sim,i,j,k,l,m,n),])

          }
        }
      }
    }
  }
}
}

#saveRDS(Results,"Simulation_Study_Breakdown_Test_Corrected.rds")
Results <- readRDS("Simulation_Study_Breakdown_Test_Corrected.rds")

library(inlabru)

# Change the factors to ordered factors to improve plotting
Results$sat_level <- factor(Results$sat_level, levels=c('low','high other', 'high target'), ordered = T)
Results$mean_attract <- factor(Results$mean_attract, levels=c('constant','linear target', 'linear aggressive'), ordered = T)
Results$correlation <- factor(Results$correlation, levels=c('negative','none', 'positive'), ordered = T)
Results$target_species <- factor(Results$target_species, levels=c('slow','equally_matched', 'aggressive'), ordered = T)
Results$sd_species <- factor(Results$sd_species, levels=c('equal','overdisp_target', 'overdisp_other'), ordered = T)
Results$sat_effects <- factor(Results$sat_effects, levels=c('target','no saturation'), ordered = T)

Results %>%
  mutate(Tested = ifelse(target_species!='aggressive' &
                           sd_species == 'overdisp_other',
                         'Experiments 2-4',
                         ifelse(target_species=='aggressive' &
                                  sd_species == 'overdisp_other',
                                'Experiment 5','Untested'))) %>%
  group_by(target_species, sd_species) %>%
  ggplot(aes(x=sd_species, y=Delta_Mean_Catch_Count, 
             group=sd_species, colour=Tested)) +
  geom_boxplot(position = position_dodge(width=0.8)) +
  scale_fill_brewer(palette = 'Pastel1') +
  geom_hline(yintercept = 0) +
  ylab('Estimated change in log mean catch count of target species (i.e. \u03b3 - \u03b2)') +
  facet_grid( ~ target_species , scales = 'free_y',
              labeller = labeller(
                target_species=c(
                  slow = 'Target species slow to reach baits',
                  equally_matched = 'Target species equally matched to reach baits',
                  aggressive = 'Target species reaches bait quickly'
                )
              )
  ) +
  xlab('Which species exhibits strongest schooling behaviour (i.e. most overdispersed catch counts)?') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.8,0.2),
        legend.title = element_blank())+#'none') +
  scale_x_discrete(labels= c('None','Target','Non-target')) +
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  geom_hline(yintercept = log(0.9), colour='black', linetype='dashed') #+ 
  #geom_hline(yintercept = log(0.8), colour='black', linetype='dotted')

# read in the plot from the sampled trajectories
plot <- readRDS('Sim_Exp_6_Trajectories_plot.rds')

multiplot(plot,
Results %>%
  mutate(Tested = factor(
    ifelse(target_species!='aggressive' &
             sd_species == 'overdisp_other',
           'red',
           ifelse(target_species=='aggressive' &
                    sd_species == 'overdisp_other',
                  'green',
                  ifelse(target_species=='aggressive' &
                           sd_species == 'overdisp_target',
                         'blue','black'))),
    levels=c('red','green','blue', 'black'),
    ordered=T),
    Tested2 = factor(
      ifelse(target_species!='aggressive' &
               sd_species == 'overdisp_other',
             'solid',
             ifelse(target_species=='aggressive' &
                      sd_species == 'overdisp_other',
                    'dotted','longdashed')),
      levels=c('solid','dotted','longdashed'),
      ordered=T),
    sd_species=factor(sd_species, levels=c('overdisp_other','equal','overdisp_target'), ordered=T)) %>%
  group_by(target_species, sd_species) %>%
  ggplot(aes(x=target_species, y=Delta_Mean_Catch_Count, 
             group=target_species, colour=Tested, linetype=Tested2)) +
  geom_boxplot(position = position_dodge(width=0.8)) +
  scale_color_discrete(guide = 'none',type=c('#F8766D', '#00BA38', '#619CFF','#000000')) +
  scale_linetype_manual(values=c('solid','dotted','longdash'),
                        labels=c('Experiments 2-4','Experiment 5','Untested')) +
  geom_hline(yintercept = 0) +
  ylab('Estimated change in log mean catch count of target species') + #(i.e. \u03b3 - \u03b2)') +
  facet_grid( ~ sd_species , scales = 'free_y',
              labeller = labeller(
                sd_species=c(
                  overdisp_other = 'A non-target species exhibits strongest \nschooling behaviour',
                  equal = 'No species exhibits schooling \nbehaviour',
                  overdisp_target = 'Target species exhibits strongest \nschooling behaviour'
                  
                )
              )
  ) +
  xlab('Which species reaches the baits fastest on average?') + guides(fill='none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.15,0.8),
        legend.key.size = unit(1,'cm'),
        legend.title = element_blank())+#'none') +
  scale_x_discrete(labels= c('Non-target','No species','Target')) +
  geom_hline(yintercept = log(0.9), colour='pink', linetype='dashed') 
)

Results %>%
  mutate(
         Tested = factor(
           ifelse(target_species!='aggressive' &
                           sd_species == 'overdisp_other',
                         'solid',
                         ifelse(target_species=='aggressive' &
                                  sd_species == 'overdisp_other',
                                'dotted','longdashed')),
         levels=c('solid','dotted','longdashed'),
         ordered=T))%>%
  group_by(target_species, sd_species) %>%
  ggplot(aes(x=sd_species, y=Delta_Mean_Catch_Count, 
             group=sd_species, colour=Tested, linetype=Tested)) +
  geom_boxplot(position = position_dodge(width=0.8)) +
  scale_fill_brewer(palette = 'Pastel1') +
  geom_hline(yintercept = 0) +
  ylab('Estimated change in log mean catch count of target species (i.e. \u03b3 - \u03b2)') +
  facet_grid( ~ target_species , scales = 'free_y',
              labeller = labeller(
                target_species=c(
                  slow = 'Target species slow to reach baits',
                  equally_matched = 'Target species equally matched to reach baits',
                  aggressive = 'Target species reaches bait quickly'
                )
              )
  ) +
  xlab('Which species exhibits strongest schooling behaviour (i.e. most overdispersed catch counts)?') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        legend.position = c(0.8,0.2),
        legend.title = element_blank())+#'none') +
  scale_x_discrete(labels= c('None','Target','Non-target')) +
  scale_color_discrete(labels=c('Experiments 2-4','Experiment 5', 'Untested')) +
  scale_linetype_discrete(labels=c('Experiments 2-4','Experiment 5','Untested')) +
  guides(color=guide_legend(override.aes=list(fill=NA))) + 
  geom_hline(yintercept = log(0.9), colour='black', linetype='dashed') #+ 
#geom_hline(yintercept = log(0.8), colour='black', linetype='dotted')


Results %>%
  filter(correlation == 'none',
         sat_effects == 'no saturation',
         target_species=='slow') %>%
  group_by(sat_level, sat_effects, mean_attract, correlation, target_species, sd_species) %>%
  mutate(Mean = median(Delta_Mean_Catch_Count),
         UCL = quantile(Delta_Mean_Catch_Count, prob=1),
         LCL = quantile(Delta_Mean_Catch_Count, prob=0)) %>%
  mutate(random_samp = c(T, rep(F, length(Delta_Mean_Catch_Count)-1))) %>%
  ggplot(aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, group=target_species)) +
  # geom_rect(data= ~.x[.x$random_samp==T,],
  #           aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, colour=target_species, group=target_species, fill = mean_attract),
  #           xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_errorbar(position = position_dodge(width=0.8)) +
  geom_point(position = position_dodge(width=0.8), size=2) +
  scale_fill_brewer(palette = 'Pastel1') +
  geom_hline(yintercept = 0) +
  #ggtitle('Model-estimated change in log mean catch counts as hook saturation increased from 85-95% to > 95%')+#,
          #subtitle = 'Values above the line indicate mean catch count increases at high levels of hook saturation, values below indicate it decreases.\nPlotted are the median, max, and minimum values seen in the 20 simulation replicates') +
  ylab('Estimated change in log mean catch count (i.e. \u03b3 - \u03b2)') +
  facet_grid( mean_attract ~ sat_level , scales = 'free_y',
              labeller = labeller(
                mean_attract = c(
                  constant = 'No change in any species\' \nabundance',
                  `linear target` = 'Target species\' abundance \nincreasing',
                  `linear aggressive` = 'Non-target species\' abundance \nincreasing '
                ),
                sat_level = c(
                  low = 'Saturation less "common"',
                  `high other` = 'Saturation "common" due to \nhigh abundance of non-target species',
                  `high target` = 'Saturation "common" due to \nhigh abundance of target species'
                )
              )) +
  xlab('Which species exhibit schooling behaviour (i.e. overdispersed catch counts)?') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none') +
  # scale_color_viridis_d(labels=c('Slow','Equally Matched','Aggressive')) +
  # scale_shape_manual(labels=c('Slow','Equally Matched','Aggressive'),
  #                    values=c('circle','triangle','square')) +
  # labs(shape = 'Target Species\' \nbehaviour',
  #      colour= 'Target Species\' \nbehaviour') +
  scale_x_discrete(labels= c('None','Target','Non-target')) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

Results %>%
  filter(correlation == 'none',
         sat_effects == 'no saturation',
         target_species=='aggressive') %>%
  group_by(sat_level, sat_effects, mean_attract, correlation, target_species, sd_species) %>%
  mutate(Mean = median(Delta_Mean_Catch_Count),
         UCL = quantile(Delta_Mean_Catch_Count, prob=1),
         LCL = quantile(Delta_Mean_Catch_Count, prob=0)) %>%
  mutate(random_samp = c(T, rep(F, length(Delta_Mean_Catch_Count)-1))) %>%
  ggplot(aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, group=target_species)) +
  # geom_rect(data= ~.x[.x$random_samp==T,],
  #           aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, colour=target_species, group=target_species, fill = mean_attract),
  #           xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_errorbar(position = position_dodge(width=0.8)) +
  geom_point(position = position_dodge(width=0.8), size=2) +
  scale_fill_brewer(palette = 'Pastel1') +
  geom_hline(yintercept = 0) +
  #ggtitle('Model-estimated change in log mean catch counts as hook saturation increased from 85-95% to > 95%')+#,
  #subtitle = 'Values above the line indicate mean catch count increases at high levels of hook saturation, values below indicate it decreases.\nPlotted are the median, max, and minimum values seen in the 20 simulation replicates') +
  ylab('Estimated change in log mean catch count (i.e. \u03b3 - \u03b2)') +
  facet_grid( mean_attract ~ sat_level , scales = 'free_y',
              labeller = labeller(
                mean_attract = c(
                  constant = 'No change in any species\' \nabundance',
                  `linear target` = 'Target species\' abundance \nincreasing',
                  `linear aggressive` = 'Non-target species\' abundance \nincreasing '
                ),
                sat_level = c(
                  low = 'Saturation less "common"',
                  `high other` = 'Saturation "common" due to \nhigh abundance of non-target species',
                  `high target` = 'Saturation "common" due to \nhigh abundance of target species'
                )
              )) +
  xlab('Which species exhibit schooling behaviour (i.e. overdispersed catch counts)?') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none') +
  # scale_color_viridis_d(labels=c('Slow','Equally Matched','Aggressive')) +
  # scale_shape_manual(labels=c('Slow','Equally Matched','Aggressive'),
  #                    values=c('circle','triangle','square')) +
  # labs(shape = 'Target Species\' \nbehaviour',
  #      colour= 'Target Species\' \nbehaviour') +
  scale_x_discrete(labels= c('None','Target','Non-target')) +
  ylim(-0.15,NA) +
  guides(color=guide_legend(override.aes=list(fill=NA)))

Results %>%
  filter(correlation == 'none',
         sat_effects == 'no saturation',
         target_species=='equally_matched') %>%
  group_by(sat_level, sat_effects, mean_attract, correlation, target_species, sd_species) %>%
  mutate(Mean = median(Delta_Mean_Catch_Count),
         UCL = quantile(Delta_Mean_Catch_Count, prob=1),
         LCL = quantile(Delta_Mean_Catch_Count, prob=0)) %>%
  mutate(random_samp = c(T, rep(F, length(Delta_Mean_Catch_Count)-1))) %>%
  ggplot(aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, group=target_species)) +
  # geom_rect(data= ~.x[.x$random_samp==T,],
  #           aes(x=sd_species, y=Mean, ymin=LCL, ymax=UCL, colour=target_species, group=target_species, fill = mean_attract),
  #           xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf,alpha = 0.3) +
  geom_errorbar(position = position_dodge(width=0.8)) +
  geom_point(position = position_dodge(width=0.8), size=2) +
  scale_fill_brewer(palette = 'Pastel1') +
  geom_hline(yintercept = 0) +
  #ggtitle('Model-estimated change in log mean catch counts as hook saturation increased from 85-95% to > 95%')+#,
  #subtitle = 'Values above the line indicate mean catch count increases at high levels of hook saturation, values below indicate it decreases.\nPlotted are the median, max, and minimum values seen in the 20 simulation replicates') +
  ylab('Estimated change in log mean catch count (i.e. \u03b3 - \u03b2)') +
  facet_grid( mean_attract ~ sat_level , scales = 'free_y',
              labeller = labeller(
                mean_attract = c(
                  constant = 'No change in any species\' \nabundance',
                  `linear target` = 'Target species\' abundance \nincreasing',
                  `linear aggressive` = 'Non-target species\' abundance \nincreasing '
                ),
                sat_level = c(
                  low = 'Saturation less "common"',
                  `high other` = 'Saturation "common" due to \nhigh abundance of non-target species',
                  `high target` = 'Saturation "common" due to \nhigh abundance of target species'
                )
              )) +
  xlab('Which species exhibit schooling behaviour (i.e. overdispersed catch counts)?') + guides(fill='none') +
  scale_fill_brewer(palette = 'Pastel1') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none') +
  # scale_color_viridis_d(labels=c('Slow','Equally Matched','Aggressive')) +
  # scale_shape_manual(labels=c('Slow','Equally Matched','Aggressive'),
  #                    values=c('circle','triangle','square')) +
  # labs(shape = 'Target Species\' \nbehaviour',
  #      colour= 'Target Species\' \nbehaviour') +
  scale_x_discrete(labels= c('None','Target','Non-target')) +
  ylim(-0.4,NA) +
  guides(color=guide_legend(override.aes=list(fill=NA)))
