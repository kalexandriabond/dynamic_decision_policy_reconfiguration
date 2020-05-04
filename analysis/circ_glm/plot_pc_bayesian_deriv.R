
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

bic_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/
pc_bic_df.csv'))

bic_df$subj_id_anon <- factor(bic_df$subj_id_anon)

pc_alt_v_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/
pc_BF_df.csv'))
pc_alt_v_null_df$subj_id_anon <- factor(pc_alt_v_null_df$subj_id_anon)


fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')

# looking at "raw" data first ... 

# plot BICs for each model for each subject 

ggplot(bic_df, aes(model_reordered, bic, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + 
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta\\sim t$')),
    parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab('BIC') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$bic) - 20 , max(bic_df$bic) + 20))

ggsave(paste0(fig_dir, 'sub_pc_models_raw_BIC.png'))



ggplot(bic_df, aes(subj_id_anon, bic, fill=model_reordered)) + geom_bar( 
  position = "dodge", stat = "identity", aes(fill=model_reordered)) + 
  theme_minimal(base_size = 16) + 

  scale_fill_discrete(name='model', labels = c(
    parse(text=TeX('$\\theta\\sim t$')),
    parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +

  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab('BIC') + xlab('subject') + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$bic) - 20 , max(bic_df$bic) + 20))

ggsave(paste0(fig_dir, 'sub_pc_models_raw_BIC.png'))


# null-adjusted bic

ggplot(bic_df[which(bic_df$null_adj_bic != 0),], aes(model_reordered, null_adj_bic, fill=subj_id_anon)) + geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab(parse(text=TeX('$BIC_{i} - BIC_{\\theta \\sim t}$'))) + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$null_adj_bic) - 20 , max(bic_df$null_adj_bic) + 20))

ggsave(paste0(fig_dir, 'sub_pc_models_nulladj_BIC.png'))


# plot mean BICs for each model 

ggplot(bic_df, aes(model_reordered, bic))  + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta\\sim t$')),
    parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab('BIC') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$bic) - 20 , max(bic_df$bic) + 20))

ggsave(paste0(fig_dir, 'mean_pc_models_raw_BIC.png'))


ggplot(bic_df[which(bic_df$null_adj_bic != 0),], aes(model_reordered, null_adj_bic)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab(parse(text=TeX('$BIC_{i} - BIC_{\\theta \\sim t}$'))) + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$null_adj_bic) - 20 , max(bic_df$null_adj_bic) + 20))

ggsave(paste0(fig_dir, 'mean_pc_models_nulladj_BIC.png'))

# ~Bayes Factors~

# subject-level data 

# evidence for alt (BF10)

ggplot(pc_alt_v_null_df, aes(comparison_reordered, BF, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t + p1 + p2 | \\theta \\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2 | \\theta\\sim t$')))) + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('PC models v. time null') + 
  ylab('BF10') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2")

ggsave(paste0(fig_dir, 'sub_pc_models_BF10.png'))


# evidence for time null (BF01)
ggplot(pc_alt_v_null_df, aes(comparison_reordered, BF_alt, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$ \\theta \\sim t | \\theta \\sim  t + p1 + p2 $')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p2 + t \\cdot p2 $')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Time null v. PC models') + 
  ylab('BF01') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") 


ggsave(paste0(fig_dir, 'sub_pc_models_BF01.png'))



# meaned plots 
# evidence for time null (BF01)
ggplot(pc_alt_v_null_df, aes(comparison_reordered, BF_alt)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$ \\theta \\sim t | \\theta \\sim  t + p1 + p2 $')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + t \\cdot p1$')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p2 + t \\cdot p2 $')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
    parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) +
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Time null v. PC models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_pc_models_BF01.png'))



# evidence for alt (BF10)
ggplot(pc_alt_v_null_df, aes(comparison_reordered, BF)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t + p1 + p2 | \\theta \\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2 | \\theta\\sim t$')),
    parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2 | \\theta\\sim t$')))) + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('PC models v. time null') + 
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_pc_models_BF10.png'))




# scale_x_discrete labels 
# t (epoch_trial)
# p1 (projection_0_shifted)
# p2 (projection_1_shifted)
# \theta (theta_radians_z)

# \\sim (tilde)
# \\cdot (*)


# BF10: model_i against time null 

# scale_x_discrete(labels = c(
#   parse(text=TeX('$\\theta \\sim  t + p1 + p2 | \\theta \\sim t$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1 | \\theta\\sim t$')),
#   parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2 | \\theta\\sim t$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2 | \\theta\\sim t$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2 | \\theta\\sim t$')))) + 


# BF01: time null against model_i

# scale_x_discrete(labels = c(
#   parse(text=TeX('$ \\theta \\sim t | \\theta \\sim  t + p1 + p2 $')),
#   parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + t \\cdot p1$')),
#   parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p2 + t \\cdot p2 $')),
#   parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
#   parse(text=TeX('$\\theta\\sim t | \\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))) + 


# models 

# 
# scale_x_discrete(labels = c
#   parse(text=TeX('$\\theta\\sim t$')),
#   parse(text=TeX('$\\theta \\sim  t + p1 + p2$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + t \\cdot p1$')),
#   parse(text=TeX('$\\theta\\sim  t + p2 + t \\cdot p2$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1 + t \\cdot p2$')),
#   parse(text=TeX('$\\theta\\sim  t + p1 + p2 + t \\cdot p1\\cdot p2$')))



