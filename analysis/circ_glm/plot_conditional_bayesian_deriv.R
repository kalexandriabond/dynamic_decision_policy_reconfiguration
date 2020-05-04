
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


fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')


cond_alt_v_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/
conditional_BF_df.csv'))

cond_alt_v_null_df$comparison_reordered <- factor(cond_alt_v_null_df$comparison, levels=
                                                    c("theta_radians_z ~ epoch_trial : theta_radians_z ~ 1",
                                                      "theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ 1",
                                                      "theta_radians_z ~ epoch_trial + p_optimal + p_optimal * epoch_trial : theta_radians_z ~ 1",
                                                      "theta_radians_z ~ epoch_trial + lambda_val + lambda_val * epoch_trial : theta_radians_z ~ 1",
                                                      "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ 1",
                                                      "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ 1"))



cond_alt_v_null_df$subj_id_anon <- factor(cond_alt_v_null_df$subj_id_anon)

# check level ordering 
print(levels(cond_alt_v_null_df$comparison_reordered))


cond_alt_v_time_null_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/
conditional_timenull_BF_df.csv'))

cond_alt_v_time_null_df$subj_id_anon <- factor(cond_alt_v_time_null_df$subj_id_anon)

cond_alt_v_time_null_df$comparison_reordered <- factor(cond_alt_v_time_null_df$comparison, levels=
                                                         c("theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ epoch_trial"))
                                                       
                                                       
                                                       
# check level ordering 
print(levels(cond_alt_v_time_null_df$comparison_reordered))



bic_df <- read.csv(paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/
conditional_bic_df.csv'))

bic_df$subj_id_anon <- factor(bic_df$subj_id_anon)

bic_df$model_reordered <- factor(bic_df$model, levels=c("abs_null_m", 
                                                        "time_null_m",
                                                        "conditional_m",
                                                        "conflict_time_m",
                                                        "vol_time_m",
                                                        "conditional_time_m",
                                                        "vol_conflict_time_m"
))


# check ordering of levels 
print(levels(bic_df$model_reordered))

# looking at "raw" data first ... 

# plot BICs for each model for each subject 

ggplot(bic_df, aes(model_reordered, bic, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + 
  theme_minimal(base_size = 16) + 

  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim t$')),
    parse(text=TeX('$\\theta \\sim t + p + \\lambda$')),
    parse(text=TeX('$\\theta \\sim t + p + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + \\lambda \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot t + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +

  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab('BIC') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$bic) - 20 , max(bic_df$bic) + 20))

ggsave(paste0(fig_dir, 'sub_cond_models_raw_BIC.png'))


# null-adjusted bic

ggplot(bic_df[which(bic_df$null_adj_bic != 0),], aes(model_reordered, null_adj_bic, fill=subj_id_anon)) + geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
  parse(text=TeX('$\\theta \\sim t$')),
  parse(text=TeX('$\\theta \\sim t + p + \\lambda$')),
  parse(text=TeX('$\\theta \\sim t + p + p \\cdot t$')),
  parse(text=TeX('$\\theta \\sim t + \\lambda + \\lambda \\cdot t$')),
  parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot t + p \\cdot t$')),
  parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +

  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab(parse(text=TeX('$BIC_{i} - BIC_{\\theta \\sim 1}$'))) + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$null_adj_bic) - 20 , max(bic_df$null_adj_bic) + 20))

ggsave(paste0(fig_dir, 'sub_cond_models_nulladj_BIC.png'))


# plot mean BICs for each model 

ggplot(bic_df, aes(model_reordered, bic))  + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim t$')),
    parse(text=TeX('$\\theta \\sim t + p + \\lambda$')),
    parse(text=TeX('$\\theta \\sim t + p + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + \\lambda \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot t + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab('BIC') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$bic) - 20 , max(bic_df$bic) + 20))

ggsave(paste0(fig_dir, 'mean_cond_models_raw_BIC.png'))


ggplot(bic_df[which(bic_df$null_adj_bic != 0),], aes(model_reordered, null_adj_bic)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95)) +
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim t$')),
    parse(text=TeX('$\\theta \\sim t + p + \\lambda$')),
    parse(text=TeX('$\\theta \\sim t + p + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + \\lambda \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot t + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('All models') + 
  ylab(parse(text=TeX('$BIC_{i} - BIC_{\\theta \\sim 1}$'))) + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") + 
  coord_cartesian(ylim=c(min(bic_df$null_adj_bic) - 20 , max(bic_df$null_adj_bic) + 20))

ggsave(paste0(fig_dir, 'mean_cond_models_nulladj_BIC.png'))

# ~Bayes Factors~

# subject-level data 

# evidence for alt (BF10)

ggplot(cond_alt_v_null_df, aes(comparison_reordered, BF, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + p + p\\cdot t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + \\lambda \\cdot t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot t  + p \\cdot t| \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t| \\theta \\sim 1$')))) +
    
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Conditional models v. abs. null') + 
  ylab('BF10') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2")
ggsave(paste0(fig_dir, 'sub_cond_models_BF10.png'))


# evidence for abs null (BF01)
ggplot(cond_alt_v_null_df, aes(comparison_reordered, BF_alt, fill=subj_id_anon)) + geom_bar( 
  position = "dodge", stat = "identity") + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim 1 | \\theta \\sim  t$')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + p + p\\cdot t $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + \\lambda \\cdot t $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p + \\lambda \\cdot t  + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +
  

  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Abs. null v. conditional models') + 
  ylab('BF01') + xlab('model') + labs(fill="subject") + 
  scale_fill_brewer(palette="Dark2") 
ggsave(paste0(fig_dir, 'sub_cond_models_BF01.png'))



# meaned plots 



# evidence for alt (BF10)
ggplot(cond_alt_v_null_df, aes(comparison_reordered, BF)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + p + p\\cdot t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + \\lambda \\cdot t | \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot t  + p \\cdot t| \\theta \\sim 1$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t| \\theta \\sim 1$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Conditional models v. absolute null') + 
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_models_BF10.png'))

# evidence for time null (BF01)
ggplot(cond_alt_v_null_df, aes(comparison_reordered, BF_alt)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim 1 | \\theta \\sim  t$')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + p + p\\cdot t $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + \\lambda \\cdot t $')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p + \\lambda \\cdot t  + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim 1 |\\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +
  
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Absolute null v. conditional models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_models_BF01.png'))




# time null 

cond_alt_v_time_null_df$comparison_reordered <- factor(cond_alt_v_time_null_df$comparison, levels=
                                                         c("theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ epoch_trial",
                                                           "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ epoch_trial"))



# meaned plots 

# evidence for alt (BF10)
ggplot(cond_alt_v_time_null_df, aes(comparison_reordered, BF)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p | \\theta \\sim t$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot t + p \\cdot t | \\theta \\sim t$')),
    parse(text=TeX('$\\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t | \\theta \\sim t$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Conditional models v. time null') + 
  ylab('BF10') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_time_null_models_BF10.png'))

# evidence for time null (BF01)
ggplot(cond_alt_v_time_null_df, aes(comparison_reordered, BF_alt)) + geom_bar(stat = "summary", fun.y = "mean") + 
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",position=position_dodge(width=0.85), width=0.1, fun.args=list(conf.int=.95))  + scale_y_log10() +  
  theme_minimal(base_size = 16) + 
  
  
  scale_x_discrete(labels = c(
    parse(text=TeX('$\\theta \\sim t | \\theta \\sim  t + \\lambda + p $')),
    parse(text=TeX('$\\theta \\sim t | \\theta \\sim  t + \\lambda + p + \\lambda \\cdot t + p \\cdot t$')),
    parse(text=TeX('$\\theta \\sim t | \\theta \\sim  t + \\lambda + p + \\lambda \\cdot p \\cdot t$')))) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept=1)  + ggtitle('Time null v. conditional models') + 
  ylab('BF01') + xlab('model') +
  scale_fill_brewer(palette="Dark2") 

ggsave(paste0(fig_dir, 'mean_cond_time_null_models_BF01.png'))



# scale_x_discrete labels 
# t (epoch_trial)
# p (p_optimal)
# \lambda (vol)
# \theta (theta_radians_z)

# \\sim (tilde)
# \\cdot (*)


# BF10: model_i against time null 

# scale_x_discrete(labels = c(
#   parse(text=TeX('$\\theta \\sim  t  | \\theta \\sim 1$')),


#   parse(text=TeX('$\\theta \\sim t +  \\lambda +  p | \\theta \\sim 1$')),

#   parse(text=TeX('$\\theta \\sim t +  \\lambda +  \\lambda \\cdot t | \\theta \\sim 1$')),
#   parse(text=TeX('$\\theta \\sim t +  p +  p \\cdot t | \\theta \\sim 1$')),

#   parse(text=TeX('$\\theta \\sim t +  \\lambda +  p  + \\lambda \\cdot t + p \\cdot t | \\theta \\sim 1$')),
#   parse(text=TeX('$\\theta \\sim t +  \\lambda +  p  + \\lambda \\cdot  p \\cdot t | \\theta \\sim 1$')),




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


# # desired factor order in terms of complexity for absolute null model comparison 
# # cond_alt_v_null_df
# c("theta_radians_z ~ epoch_trial : theta_radians_z ~ 1",
#   "theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ 1",
#   "theta_radians_z ~ epoch_trial + p_optimal + p_optimal * epoch_trial : theta_radians_z ~ 1",
#   "theta_radians_z ~ epoch_trial + lambda_val + lambda_val * epoch_trial : theta_radians_z ~ 1",
#   "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ 1",
#   "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ 1")
#   
# # desired factor order in terms of complexity for time null model comparison 
# # cond_alt_v_time_null_df
# c("theta_radians_z ~ epoch_trial + lambda_val + p_optimal : theta_radians_z ~ epoch_trial",
#   "theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val * epoch_trial + p_optimal * epoch_trial : theta_radians_z ~ epoch_trial",
#   "theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal * lambda_val * epoch_trial : theta_radians_z ~ epoch_trial")
#   
#   
# do(abs_null_m = circGLM(theta_radians_z ~ 1, data=., Q = Q, burnin = burnin, thin = thin), 
#    
#    time_null_m = circGLM(theta_radians_z ~ epoch_trial, data=., Q = Q, burnin = burnin, thin = thin),
#    
#    conditional_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal, data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    conditional_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + p_optimal + lambda_val*epoch_trial + p_optimal*epoch_trial, 
#                                 data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    vol_time_m = circGLM(theta_radians_z ~ epoch_trial + lambda_val + lambda_val*epoch_trial, 
#                         data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    
#    conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + p_optimal*epoch_trial, 
#                              data =., Q = Q, burnin = burnin, thin = thin), 
#    
#    vol_conflict_time_m = circGLM(theta_radians_z ~ epoch_trial + p_optimal + lambda_val + p_optimal*lambda_val*epoch_trial, 
#                                  data =., Q = Q, burnin = burnin, thin = thin))
# 

  