# written by KAMB circa 8/2020

rm(list=ls())

library(dplyr)
library(broom)
library(tidyverse)
library(circglmbayes)
library(reshape2)
library(boot)
library(latex2exp)


if (Sys.info()['sysname'] == 'Darwin'){
  home = '/Users/67981492/'
} else if (Sys.info()['sysname'] == 'Linux'){
  home = '/data/'
}


# here, BF = model_n / null (BF10) & BF_alt = null / model_n (BF01)

# comparing to the time null model (BF01: time / model_n)
pc_alt_v_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_BF_df.csv'))
pc_alt_v_null_df$subj_id_anon <- factor(pc_alt_v_null_df$subj_id_anon)


fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')

# model naming to estimation reference 

# pc0 = timing component; pc1 = mass component

model_comparison_order <- c("theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted : theta_radians_z ~ epoch_trial",
                 "theta_radians_z ~ epoch_trial + projection_0_shifted + epoch_trial * projection_0_shifted : theta_radians_z ~ epoch_trial",
                  "theta_radians_z ~ epoch_trial + projection_1_shifted + epoch_trial * projection_1_shifted : theta_radians_z ~ epoch_trial",
                 "theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + epoch_trial * projection_0_shifted + epoch_trial * projection_1_shifted : theta_radians_z ~ epoch_trial",
                 "theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + projection_0_shifted * projection_1_shifted * epoch_trial : theta_radians_z ~ epoch_trial")
                 

# meaned plots 

# evidence for alt (BF10)
ggplot(pc_alt_v_null_df, aes(comparison, BF)) + geom_bar(stat = "summary", fun = mean) + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 20) + 
  scale_x_discrete(limits = model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  geom_hline(yintercept=1)  + ggtitle('PC models v. time null') + 
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_pc_time_null_models_BF10_raw_labels.png'))


# evidence for time null (BF01)
ggplot(pc_alt_v_null_df, aes(comparison, BF_alt)) + geom_bar(stat = "summary", fun = mean) + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 20) + 
  scale_x_discrete(limits = model_comparison_order) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  geom_hline(yintercept=1)  + ggtitle('Time null v. PC models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_pc_time_null_models_BF01_raw_labels.png'))






