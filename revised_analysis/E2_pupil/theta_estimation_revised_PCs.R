# written by KAMB circa 8/2020


# setwd('~/Dropbox/loki_0.5/analysis/circ_glm/')


# note that this is updated to include additional PCs (first 4)

rm(list=ls())

library(dplyr)
library(broom)
library(tidyverse)
library(circglmbayes)
library(reshape2)
library(boot)


# if (Sys.info()['sysname'] == 'Darwin'){
#   home = '/Users/67981492/'
# } else if (Sys.info()['sysname'] == 'Linux'){
#   home = '/data/'
# }


av_manifold_df = read.csv('Documents/elife_revisions_loki/pupil_manifold_df.csv', header = TRUE)
head(av_manifold_df)

av_manifold_df_sorted <- av_manifold_df[order(av_manifold_df$subj_id, av_manifold_df$condition, av_manifold_df$trial),]
av_manifold_df_sorted_subset <- subset(av_manifold_df_sorted, shifted_epoch_trial <= 3 & shifted_epoch_trial >= 0)


# av_manifold_df_sorted_subset$PC0 = pc_ls_df_sorted_subset$PC0
# av_manifold_df_sorted_subset$PC1 = pc_ls_df_sorted_subset$PC1
# # av_manifold_df_sorted_subset$projection_2 = pc_ls_df_sorted_subset$projection_2

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

ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, PC0)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, PC1)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)


ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, PC2)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, PC3)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)



ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), PC0)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), PC1)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), PC2)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(lambda_val), PC3)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" ) + 
  facet_wrap(~ subj_id, ncol=2)


ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), PC0)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), PC1)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), PC2)) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2)

ggplot(av_manifold_df_sorted_subset, aes(factor(p_optimal), PC3)) + 
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
ggplot(av_manifold_df_sorted_subset, aes(lambda_val, epoch_number, color=factor(lambda_val))) + 
  stat_summary(fun.data = mean_cl_boot, geom = "pointrange" )  + 
  facet_wrap(~ subj_id, ncol=2) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 


# hues for interactions
ggplot(av_manifold_df_sorted_subset, aes(epoch_trial, color=factor(lambda_val))) + 
  geom_histogram() + 
  facet_wrap(~ subj_id, ncol=2) 
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 



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



# models testing hypotheses regarding principal component effects on theta

pc_models = av_manifold_df_sorted_subset_clean %>% group_by(subj_id) %>%
  
  do(null_m = circGLM(theta_radians_z ~ epoch_trial, data=.,  Q = Q, burnin = burnin, thin = thin), 
    
     
     null_pc_interaction_m = circGLM(theta_radians_z ~ epoch_trial + PC0 + PC1 + PC2 + PC3, data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     
     pc0123time_interaction_m = circGLM(theta_radians_z ~ epoch_trial + PC0 + PC1 + PC2 + PC3 + epoch_trial*PC0 + 
                               epoch_trial*PC1  + epoch_trial*PC2 + epoch_trial*PC3, data =.,  Q = Q, burnin = burnin, thin = thin), 
     
     
     pc0_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + PC0 + epoch_trial*PC0, data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     pc1_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + PC1 + epoch_trial*PC1, data =.,  Q = Q, burnin = burnin, thin = thin),
     
     pc2_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + PC2 + epoch_trial*PC2, data =., Q = Q, burnin = burnin, thin = thin), 
     
     
     pc3_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + PC3 + epoch_trial*PC3, data =.,  Q = Q, burnin = burnin, thin = thin),
    
     pc01time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + PC0 + PC1 + PC2 + PC3 + PC0*epoch_trial + PC1*epoch_trial + PC2*epoch_trial + PC3*epoch_trial, data =.,  Q = Q, burnin = burnin, thin = thin)
     
     
     )

print(paste('done with pc models. saving workspace image ... '))


save.image(file='/Users/i_67981492/Documents/elife_revisions_loki/circ_glms_pupil_deriv_PCs.RData')
