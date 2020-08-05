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

all_models_posterior_prob_df = read.csv(file=paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/cond_postprob_each_m.csv'))
all_models_posterior_prob_df$subj_id <- as.factor(all_models_posterior_prob_df$subj_id)

# model naming to estimation reference 

# do(abs_null_m = circGLM(theta_radians_z ~ 1, data=., Q = Q, burnin = burnin, thin = thin), 
#    
#    time_null_m = circGLM(theta_radians_z ~ epoch_trial, data=., Q = Q, burnin = burnin, thin = thin),
#    
#    conditional_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal, data =., Q = Q, burnin = burnin, thin = thin), 

#    
#    vol_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + lambda_val*epoch_trial, 
#                         data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + p_optimal*epoch_trial, 
#                              data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    conditional_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val*epoch_trial + p_optimal*epoch_trial, 
#                                 data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    vol_conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal*lambda_val*epoch_trial, 
#                                  data =., Q = Q, burnin = burnin, thin = thin))

model_order <- c("abs_null_m", "time_null_m", "conditional_m",  "conflict_time_m", "vol_time_m", "conditional_time_m", "vol_conflict_time_m")


ggplot(all_models_posterior_prob_df, aes(model_reordered, posterior_prob, fill=subj_id)) + geom_bar( position = "dodge", stat = "identity") + 
ylab('p(M|D)') +   geom_hline(yintercept=c(0.50, 0.75, 0.95, 0.99), color='darkgray', linetype="dashed")  + theme_minimal(base_size = 20) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
scale_x_discrete(limits = model_order) +
scale_fill_brewer(palette="Dark2") +  labs(fill="subject") 

ggsave(paste0(fig_dir, 'sub_conditional_models_postprob.png'))


ggplot(all_models_posterior_prob_df, aes(model_reordered, posterior_prob)) + 
geom_bar(stat = "summary", fun = mean) +
stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
ylab('p(M|D)') +   geom_hline(yintercept=c(0.50, 0.75, 0.95, 0.99), color='darkgray', linetype="dashed")  + theme_minimal(base_size = 20) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
scale_x_discrete(limits = model_order) + 
scale_fill_brewer(palette="Dark2") +  labs(fill="subject") 

ggsave(paste0(fig_dir, 'mean_conditional_models_postprob.png'))

