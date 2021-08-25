rm(list=ls())


library(boot)
library(latex2exp)
library(tidyverse)

if (Sys.info()['sysname'] == 'Darwin'){
  home = '/Users/67981492/'
} else if (Sys.info()['sysname'] == 'Linux'){
  home = '/data/'
}

fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')

all_models_posterior_prob_df = read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_postprob_each_m.csv'))
all_models_posterior_prob_df$subj_id = as.factor(all_models_posterior_prob_df$subj_id)

# model naming to estimation reference 

# pc0 = timing component; pc1 = mass component

# do(null_m = circGLM(theta_radians_z ~ epoch_trial, data=.,  Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    null_pc_interaction_m = circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted, data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    
#    pc0time_pc1time_interaction_m = circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + epoch_trial*projection_0_shifted + 
#                                              epoch_trial*projection_1_shifted, data =.,  Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    pc0_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + epoch_trial*projection_0_shifted, data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    pc1_time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_1_shifted + epoch_trial*projection_1_shifted, data =.,  Q = Q, burnin = burnin, thin = thin),
#    
#    pc01time_interaction_m =  circGLM(theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + projection_0_shifted*projection_1_shifted*epoch_trial, data =.,  Q = Q, burnin = burnin, thin = thin)
#    
#    
# )


model_order <- c("null_m", "null_pc_interaction_m", "pc0_time_interaction_m", "pc1_time_interaction_m",
               "pc0time_pc1time_interaction_m", "pc01time_interaction_m")

ggplot(all_models_posterior_prob_df, aes(model, posterior_prob, fill=subj_id)) + geom_bar( position = "dodge", stat = "identity") + 
ylab('p(M|D)') +   geom_hline(yintercept=c(0.50, 0.75, 0.95, 0.99), color='darkgray', linetype="dashed")  + theme_minimal(base_size = 20) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
scale_x_discrete(limits = model_order) +
scale_fill_brewer(palette="Dark2") +  labs(fill="subject") 

ggsave(paste0(fig_dir, 'sub_pc_models_postprob.png'))


ggplot(all_models_posterior_prob_df, aes(model, posterior_prob)) + geom_bar(stat = "summary", fun = mean) + 
stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
ylab('p(M|D)') +   geom_hline(yintercept=c(0.50, 0.75, 0.95, 0.99), color='darkgray', linetype="dashed")  + theme_minimal(base_size = 20) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_x_discrete(limits = model_order) +
scale_fill_brewer(palette="Dark2") +  labs(fill="subject") 

ggsave(paste0(fig_dir, 'mean_pc_models_postprob.png'))