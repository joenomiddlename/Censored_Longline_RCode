# Rscript for viewing how the mean catch count changes with proportion of bait removed
# View trajectories for simulated data
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)
seed <- 25042021 # date
set.seed(seed)
nspecies <- 3
bite_funs <-cbind(c('exp_decay','mixed', 'constant'),
                  c('constant','mixed','exp_decay')
)
soak_time <- 5
n_hooks <- 800
nstation <- 6
nyears <- 100
saturation_level <- c(2)
mean_attract <- c('constant')
cov_matrices <- list( overdisp=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,0,0,0,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2))
                      )
n_sim <- 1
hook_sat_level <- 0.85 # true proportion at which saturation effects begin
cprop=0.85 # assumed proportion at which saturation effects begin
#cprop2=0.95 # assumed proportion "" "" second model
upper_bound_quantiles <- c(1)


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

Results <- vector('list',length = 3)
names(Results) <- c('Schooling Aggressive without Schooling Competitor',
                    'Non-Schooling Slow with Schooling Competitor',
                    'Non-schooling Aggressive with Schooling Competitor')

j=1
      for(k in 1:dim(bite_funs)[2])
      {
        i=1
        # NOTE THAT WE NOW SCALE THE MEAN OF THE TARGET SPECIES BY THE EXPONENTIAL OF THE
        # LOG-NORMAL VARIANCE TO ENSURE THE SAME MEAN FROM THE PREVIOUS EXPERIMENT
        # DOESN'T AFFECT relative abundance!!
        if(mean_attract[j] == 'constant')
        {
          mean_bite_gen =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (rep(100, 6)/exp(cov_matrices[[1]][3,3]/2)))
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
                  (c(120, 140, 160, 180, 200, 220)-100)/exp(cov_matrices[[1]][3,3]/2))
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  (c(120, 140, 160, 180, 200, 220)-100))
        }
        # for(l in 1:dim(bite_funs)[2])
        # {
        saturation_effect <- c(0.8, 0.8, 0.8)
        # sample the number of each species that WOULD bite at each station for each year if hooks were available
        nbite <- data.frame(bites = rep(0, times=nspecies*nstation*nyears),
                            attracted = rep(0, times=nspecies*nstation*nyears),
                            species=rep(1:3, each=nstation*nyears),
                            station=rep(1:nstation, times=nspecies*nyears),
                            year=rep(rep(1:nyears, each=nstation), times=nspecies))
        nbite$attracted <- rpois(dim(nbite)[1],
                                 lambda = as.numeric(mean_bite_gen[cbind(nbite$station,nbite$species)])*
                                   exp(as.numeric(rmvn(n=nstation*nyears, mu = rep(0,nspecies), V=cov_matrices[[1]])[cbind(rep(1:(nstation*nyears),times=nspecies),nbite$species)])))
        #exp(rnorm(dim(nbite)[1], mean = 0, sd=sd_log_bite[nbite$species])))
        
        for(i2 in 1:nstation)
        {
          for(j2 in 1:nyears)
          {
            bite_time_1 <- bite_samp(bite_funs[1,k],sum(nbite$attracted[nbite$species==1 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_1)>soak_time)
            {
              bite_time_1[bite_time_1>soak_time] <-
                bite_samp(bite_funs[1,k],sum(bite_time_1>soak_time))
            }
            # repeat for species 2 from uniform distribution
            bite_time_2 <- bite_samp(bite_funs[2,k],sum(nbite$attracted[nbite$species==2 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_2)>soak_time)
            {
              bite_time_2[bite_time_2>soak_time] <-
                bite_samp(bite_funs[2,k],sum(bite_time_2>soak_time))
            }
            
            # repeat for species 3 from uniform distribution
            bite_time_3 <- bite_samp(bite_funs[3,k],sum(nbite$attracted[nbite$species==3 &
                                                                        nbite$station==i2 &
                                                                        nbite$year==j2]))
            # truncate them to 0-5 interval
            while(max(bite_time_3)>soak_time)
            {
              bite_time_3[bite_time_3>soak_time] <-
                bite_samp(bite_funs[3,k],sum(bite_time_3>soak_time))
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
        
        if(k==1)
        {
          # The competitor species is overdispersed/schooling and fast
          # The target species is slow and non-schooling
          # This matches the sim exp 2-4 settings. 
          # Store the catch counts of both the competitor and the target
          Results[[1]] <- nbite[nbite$species==1,]
          Results[[2]] <- nbite[nbite$species==3,]
        }
        if(k==2)
        {
          # The competitor species is overdispersed/schooling but slow
          # The target species is fast and non-schooling
          # This matches the sim exp 5 setting. 
          # Store the catch counts of the target
          Results[[3]] <- nbite[nbite$species==3,]
        }
        
      }

# Plot the trajectories using the same mgcv model as in case study

mod <-
  mgcv::gam(
    formula= bites ~ -1 + s(prop_sat) + factor(year),
    family='nb',
    data=Results[[1]]
  )
mod2 <-
  mgcv::gam(
    formula= bites ~ -1 + s(prop_sat) + factor(year),
    family='nb',
    data=Results[[2]]
  )
mod3 <-
  mgcv::gam(
    formula= bites ~ -1 + s(prop_sat) + factor(year),
    family='nb',
    data=Results[[3]]
  )

pred_df <-
  expand.grid(prop_sat=seq(from=0.6,to=1, by=0.005),
              year=1:6)

pred_mod <-
  predict.gam(
    mod, 
    newdata = pred_df,
    se.fit = T,
    type = 'terms',
    terms='s(prop_sat)'
  )
pred_mod2 <-
  predict.gam(
    mod2, 
    newdata = pred_df,
    se.fit = T,
    type = 'terms',
    terms='s(prop_sat)'
  )
pred_mod3 <-
  predict.gam(
    mod3, 
    newdata = pred_df,
    se.fit = T,
    type = 'terms',
    terms='s(prop_sat)'
  )

pred_df_all <-
  rbind(cbind(pred_df,Species=1, Mean=0, LCL=0, UCL=0), 
        cbind(pred_df,Species=2, Mean=0, LCL=0, UCL=0), 
        cbind(pred_df,Species=3, Mean=0, LCL=0, UCL=0))

pred_df_all[pred_df_all$Species==1,
  c('Mean','LCL','UCL')
] <-
  c(pred_mod$fit,
    pred_mod$fit - 1.96*pred_mod$se.fit,
    pred_mod$fit + 1.96*pred_mod$se.fit
  )
pred_df_all[pred_df_all$Species==2,
         c('Mean','LCL','UCL')
] <-
  c(pred_mod2$fit,
    pred_mod2$fit - 1.96*pred_mod2$se.fit,
    pred_mod2$fit + 1.96*pred_mod2$se.fit
  )
pred_df_all[pred_df_all$Species==3,
            c('Mean','LCL','UCL')
  ] <-
  c(pred_mod3$fit,
    pred_mod3$fit - 1.96*pred_mod3$se.fit,
    pred_mod3$fit + 1.96*pred_mod3$se.fit
  )

# Plot
pred_df_all %>%
  mutate(Species=factor(
           ifelse(Species==1,
                        'Schooling \nCompetitor species of experiments 2-4',
                        ifelse(Species==2,
                        'Non-schooling with schooling competitor \nTarget species of experiments 2-4',
                        'Non-schooling quick to reach baits \nTarget species of experiment 5')),
           levels=c('Non-schooling with schooling competitor \nTarget species of experiments 2-4',
                    'Non-schooling quick to reach baits \nTarget species of experiment 5',
                    'Schooling \nCompetitor species of experiments 2-4'),
           ordered = T),
         Species=factor(Species)) %>%
  group_by(Species) %>%
  mutate(UCL = UCL-mean(Mean),
         LCL = LCL-mean(Mean),
         Mean = Mean-mean(Mean),) %>%
  ggplot(aes(x=prop_sat, y=Mean, ymax=UCL, ymin=LCL, 
             colour=Species, fill=Species)) +
  geom_rect(xmin=0.8, xmax=0.95, ymin=-1.5, ymax=1.5, colour='grey',alpha=0.2, fill='grey') +
  geom_rect(xmin=0.95, xmax=1, ymin=-1.5, ymax=1.5, colour='dark grey',alpha=0.2, fill='dark grey') +
  geom_text(data=data.frame(x=c(0.875,0.7, 0.975, hook_sat_level+0.025), y=c(-1,-1,-1,0.8), label=c('\u03B2','\u03B1',  '\u03B3','p*'), size=4),aes(x=x,y=y,label=label,size=size), inherit.aes = F) +
  geom_line() + 
  geom_ribbon(alpha=0.3) +
  #geom_vline(xintercept=0.875) +
  geom_vline(xintercept=hook_sat_level) +
  facet_grid(~Species) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  scale_color_discrete() +
  scale_fill_discrete() +                                                                
  xlab('Proportion of baits removed') +
  ylab('Model-estimated residual effect of hook competition on log mean catch count')

plot <-
pred_df_all %>%
  mutate(Species=factor(
    ifelse(Species==1,
           'Schooling \nCompetitor species of experiments 2-4',
           ifelse(Species==2,
                  'Non-schooling with schooling competitor \nTarget species of experiments 2-4',
                  'Non-schooling quick to reach baits \nTarget species of experiment 5')),
    levels=c('Non-schooling with schooling competitor \nTarget species of experiments 2-4',
             'Non-schooling quick to reach baits \nTarget species of experiment 5',
             'Schooling \nCompetitor species of experiments 2-4'),
    ordered = T),
    Species=factor(Species)) %>%
  group_by(Species) %>%
  mutate(UCL = UCL-mean(Mean),
         LCL = LCL-mean(Mean),
         Mean = Mean-mean(Mean),) %>%
  ggplot(aes(x=prop_sat, y=Mean, ymax=UCL, ymin=LCL, 
             colour=Species, fill=Species)) +
  geom_rect(xmin=0.8, xmax=0.95, ymin=-1.5, ymax=1.5, colour='grey',alpha=0.2, fill='grey') +
  geom_rect(xmin=0.95, xmax=1, ymin=-1.5, ymax=1.5, colour='dark grey',alpha=0.2, fill='dark grey') +
  geom_text(data=data.frame(x=c(hook_sat_level+0.025), y=c(0.8), label=c('p*'), size=4),aes(x=x,y=y,label=label,size=size), inherit.aes = F) +
  geom_line() + 
  geom_ribbon(alpha=0.3) +
  #geom_vline(xintercept=0.875) +
  geom_vline(xintercept=hook_sat_level) +
  facet_grid(~Species) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  scale_color_discrete() +
  scale_fill_discrete() +                                                                
  xlab('Proportion of baits removed') +
  ylab('Model-estimated residual effect of hook competition on log mean catch count')
plot

# save the plot as an rds file to merge with Sim exp 6 plots
saveRDS(plot,'Sim_Exp_6_Trajectories_plot.rds')

# Play around with the 'hook_sat_level' value (which equals p*)
# Rerun the code and notice how the point at which the curve begins to decrease
# is close to p*.