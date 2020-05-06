# setwd('~/Dropbox/loki_0.5/analysis/circ_glm/')

rm(list=ls())

library(dplyr)
library(broom)
library(tidyverse)
library(circglmbayes)
library(reshape2)
library(boot)


if (Sys.info()['sysname'] == 'Darwin'){
  home = '/Users/67981492/'
} else if (Sys.info()['sysname'] == 'Linux'){
  home = '/data/'
}


av_manifold_df = read.csv(paste0(home, "Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/aggregated_data/av_est.csv"), header = TRUE)

paste('min epoch trial ', min(av_manifold_df$shifted_epoch_trial), 'max epoch trial ',max(av_manifold_df$shifted_epoch_trial))

head(av_manifold_df)
str(av_manifold_df)


# plotting the predictors against the response
ggplot(av_manifold_df, aes(shifted_epoch_trial, theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) 


ggplot(av_manifold_df, aes(factor(lambda_val), theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  


ggplot(av_manifold_df, aes(factor(p_optimal), theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  


# hues for interactions
ggplot(av_manifold_df, aes(shifted_epoch_trial, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) 


ggplot(av_manifold_df, aes(shifted_epoch_trial, theta_radians_z, color=factor(p_optimal))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) 


ggplot(av_manifold_df, aes(p_optimal, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  



# plotting the predictors against the response by sub
ggplot(av_manifold_df, aes(shifted_epoch_trial, theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  + 
  facet_wrap(~ subj_id, ncol=4) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) 


ggplot(av_manifold_df, aes(lambda_val, theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  +
facet_wrap(~ subj_id, ncol=4) 
  

ggplot(av_manifold_df, aes(p_optimal, theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  stat_summary(fun.data = mean_cl_boot, geom = "line")  +
  facet_wrap(~ subj_id, ncol=4) 
    


# Q = 10000 * 3
# burnin = 5000
# thin = 1

Q = 10000 
burnin = 2000
thin = 1


av_manifold_df_subset = subset(av_manifold_df, epoch_trial <= 8 & epoch_trial >= 0)

av_manifold_df_clean <- av_manifold_df_subset[complete.cases(av_manifold_df_subset), ]


av_manifold_df_clean$epoch_trial <- as.factor(av_manifold_df_clean$epoch_trial)


# model testing the decay of the angular distribution to uniformity 

start_time <- Sys.time()

angular_decay_m = circGLM(theta_radians_z ~ epoch_trial, data=av_manifold_df_clean, Q = Q, burnin = burnin, thin = thin)


end_time <- Sys.time()

total_time <- end_time - start_time 

paste0('time taken to estimate models in m: ', total_time / 60 )


save.image("~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/circ_glm/angular_decay.RData")

