
rm(list=ls())

library(dplyr)
library(broom)
library(tidyverse)
library(circglmbayes)
library(reshape2)
library(boot)
library(latex2exp)
library(scales)


if (Sys.info()['sysname'] == 'Darwin'){
  home = '/Users/67981492/'
} else if (Sys.info()['sysname'] == 'Linux'){
  home = '/data/'
}


fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')


# here, BF = model_n / null (BF10) & BF_alt = null / model_n (BF01)

# comparing to the absolute null model (BF01: intercept / model_n)
cond_alt_v_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/conditional_BF_df.csv'))
cond_alt_v_null_df$subj_id_anon <- factor(cond_alt_v_null_df$subj_id_anon)

  
# comparing to the time null model (BF01: time / model_n)
cond_alt_v_time_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/conditional_timenull_BF_df.csv'))
cond_alt_v_time_null_df$subj_id_anon <- factor(cond_alt_v_time_null_df$subj_id_anon)

alt_v_abs_null_model_comparison_order <- c("theta_radians_z ~ epoch_trial : theta_radians_z ~ 1",
                            "theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ 1",
                            "theta_radians_z ~ epoch_trial + p_optimal + p_optimal * epoch_trial : theta_radians_z ~ 1",
                            "theta_radians_z ~ epoch_trial + lambda_val + lambda_val * epoch_trial : theta_radians_z ~ 1",
                            "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ 1",
                            "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ 1")


# meaned plots 

# evidence for alt (BF10)
ggplot(cond_alt_v_null_df, aes(comparison, BF)) + geom_bar(stat = "summary", fun = mean) + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 20) +
  scale_x_discrete(limits = alt_v_abs_null_model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  geom_hline(yintercept=1)  + ggtitle('Conditional models v. absolute null') + 
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_abs_null_models_BF10_raw_labels.png'))

# evidence for absolute null (BF01)
ggplot(cond_alt_v_null_df, aes(comparison, BF_alt)) + geom_bar(stat = "summary", fun = mean) + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 20) + 
  scale_x_discrete(limits = alt_v_abs_null_model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  geom_hline(yintercept=1)  + ggtitle('Absolute null v. conditional models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_abs_null_models_BF01_raw_labels.png'))




# time null 

alt_v_time_null_model_comparison_order <- c("theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ epoch_trial")



# meaned plots 

# evidence for alt (BF10)
ggplot(cond_alt_v_time_null_df, aes(comparison, BF)) + geom_bar(stat = "summary", fun = mean) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +
  theme_minimal(base_size = 20) +
  scale_x_discrete(limits = alt_v_time_null_model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  geom_hline(yintercept=1)  + ggtitle('Conditional models v. time null') +
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2")

ggsave(paste0(fig_dir, 'mean_cond_time_null_models_BF10_raw_labels.png'))

# evidence for time null (BF01)
ggplot(cond_alt_v_time_null_df, aes(comparison, BF_alt)) + 
geom_bar(stat = "summary", fun = mean ) + # error summarized in next fn with mean_cl_boot
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 20) + 
  scale_x_discrete(limits = alt_v_time_null_model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  geom_hline(yintercept=1)  + ggtitle('Time null v. conditional models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_time_null_models_BF01_raw_labels.png'))


