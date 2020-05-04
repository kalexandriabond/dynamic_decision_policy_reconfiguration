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


av_manifold_df = read.csv(paste0(home, "Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv"), header = TRUE)
head(av_manifold_df)

pc_ls_df = read.csv(paste0(home, "Dropbox/loki_0.5/analysis/aggregated_data/pc_ls_reg_est_df.csv"), header = TRUE)
head(pc_ls_df)

pc_ls_df_sorted <- pc_ls_df[order(pc_ls_df$subj_id, pc_ls_df$condition, pc_ls_df$trial),]

pc_ls_df_sorted$projection_0_shifted = lag(pc_ls_df_sorted$projection_0, n = 1)
pc_ls_df_sorted$projection_1_shifted = lag(pc_ls_df_sorted$projection_1, n = 1)


pc_ls_df_sorted_subset <- subset(pc_ls_df_sorted, shifted_epoch_trial <= 3 & shifted_epoch_trial >= 0)

av_manifold_df_sorted <- av_manifold_df[order(av_manifold_df$subj_id, av_manifold_df$condition, av_manifold_df$trial),]
av_manifold_df_sorted_subset <- subset(av_manifold_df_sorted, shifted_epoch_trial <= 3 & shifted_epoch_trial >= 0)


av_manifold_df_sorted_subset$projection_0_shifted = pc_ls_df_sorted_subset$projection_0_shifted
av_manifold_df_sorted_subset$projection_1_shifted = pc_ls_df_sorted_subset$projection_1_shifted
# av_manifold_df_sorted_subset$projection_2 = pc_ls_df_sorted_subset$projection_2

head(av_manifold_df_sorted_subset) 

str(av_manifold_df_sorted_subset)





# plotting the predictors against the response
ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange")  + 
  facet_wrap(~ subj_id, ncol=2)  + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)


ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), theta_radians_z)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)


# hues for interactions
ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z, color=factor(p_optimal))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


ggplot(av_manifold_df_sorted_subset, aes(p_optimal, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)



# pcs 

ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, projection_0_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, projection_1_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)


ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), projection_0_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), projection_1_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), projection_0_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), projection_1_shifted)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)

# # theta z by pc hues over epoch_trial 

# d + aes(colour = factor(vs)) + stat_summary(fun.y = mean, geom="line")

# ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z)) +
#   stat_summary(fun.data = mean_cl_boot,
#                geom = "pointrange") +
#   facet_wrap(~ subj_id, ncol=2)
# # 
# 
# ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_z_radians)) + 
#   stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
#   facet_wrap(~ subj_id, ncol=2)

# hues for interactions
ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, theta_radians_z, color=factor(p_optimal))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


ggplot(av_manifold_df_sorted_subset, aes(p_optimal, theta_radians_z, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)


Q = 20000
burnin = 10000
thin = 2


av_manifold_df_sorted_subset_clean <- av_manifold_df_sorted_subset[complete.cases(av_manifold_df_sorted_subset), ]

# fit model set

av_manifold_df_sorted_subset_clean$epoch_trial <- as.factor(av_manifold_df_sorted_subset_clean$epoch_trial)


# models testing hypotheses regarding conditional effects on theta


conditional_models = av_manifold_df_sorted_subset_clean %>% group_by(subj_id) %>%
  
  
  do(abs_null_m = circGLM(theta_radians_z ~ 1, data=., Q = Q, burnin = burnin, thin = thin), 
     
     time_null_m = circGLM(theta_radians_z ~ epoch_trial, data=., Q = Q, burnin = burnin, thin = thin),
     
     conditional_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal, data =., Q = Q, burnin = burnin, thin = thin), 
     
     conditional_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val*epoch_trial + p_optimal*epoch_trial, 
                                  data =., Q = Q, burnin = burnin, thin = thin), 
     
     vol_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + lambda_val*epoch_trial, 
                                  data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + p_optimal*epoch_trial, 
                          data =., Q = Q, burnin = burnin, thin = thin), 
     
     vol_conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal*lambda_val*epoch_trial, 
                                   data =., Q = Q, burnin = burnin, thin = thin))


# models testing hypotheses regarding principal component effects on theta

print(paste('done with conditional models. moving on to pc models ... '))

pc_models = av_manifold_df_sorted_subset_clean %>% group_by(subj_id) %>%
  
  do(null_m = circGLM(theta_radians_z ~ epoch_trial, data=.,  Q = Q, burnin = burnin, thin = thin), 
    
     
     null_pc_interaction_m = circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted, data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     
     pc0time_pc1time_interaction_m = circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + epoch_trial*projection_0_shifted + 
                               epoch_trial*projection_1_shifted, data =.,  Q = Q, burnin = burnin, thin = thin), 
     
     
     pc0_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + epoch_trial*projection_0_shifted, data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     pc1_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_1_shifted + epoch_trial*projection_1_shifted, data =.,  Q = Q, burnin = burnin, thin = thin),
    
     pc01time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + projection_0_shifted*projection_1_shifted*epoch_trial, data =.,  Q = Q, burnin = burnin, thin = thin)
     
     
     )

print(paste('done with pc models. saving workspace image ... '))


save.image(file=paste0(home, 'Dropbox/loki_0.5/analysis/circ_glm/circ_glms.RData'))
