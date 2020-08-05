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

load(paste0(home, 'Dropbox/loki_0.5/analysis/circ_glm/circ_glms.RData'))



# coef(models$intercept_m[[1]]) # to access a model 

# predicted = predict(pc_models$m[[1]]) # get predictions of winning model
# predfun <- predict_function.circGLM(m) # get the prediction function 

# plot_meancompare.circGLM for categorical vars 

# plot(m, type = "predict")

# get bayes factor table 
# this only includes BF for parameter hypothesis tests within a model
# does not allow arbitration between models
# BF.circGLM(m)
# BF.circGLM(interaction_m)
# BF.circGLM(intercept_m)

# # information criteria comparison of models
# # does not include BIC, but does include a Bayesian AIC 
# # which is not BIC, but Bayesian estimates of the parameters are used
# # in log-likelihood estimate during calc. of AIC

# ic_comparison = IC_compare.circGLM(m, interaction_m, intercept_m)
# ic_comparison


# so calculate the BIC 
# note that how to consider n is still an open question 
# not sure if absolute n_observations or n_independent_observations
# complicated by fitting one model for each subject, which would mean that n_independent_obs = 1
# there are papers on calculating BIC for mixed effects models but none for single-sub designs 

# books on time series analysis for other disciplines consider each obs. of X as contributing to n
# time series implies that the samples are not independent

calc_BIC <- function(model) {
  
  # calculate the BIC for a model 
  
  # the BIC is calculated for each iteration to obtain a distribution of BICS over which 
  # an empirical CI can be calculated (and can be used to calculate distributions of derivative metrics,
  # like the BF)
  
  # ll_th_estpars is the log likelihood of the model given the parameters estimated 
  # ll_th_curpars is the log likelhood of the model given the parameters estimated for mcmc iteration 
  # n_par is the number of parameters estimated for the model 
  # nrow(m$dataX) gives the number of observations used to estimate the model 
  
  # call signature: model = calc_bic(model)
  
  
  
  model$BIC = - 2 * model$ll_th_estpars + (model$n_par * log(nrow(model$data_X))) 
  model$BIC_it = -2 * model$ll_th_curpars + (model$n_par * log(nrow(model$data_X))) 
  
  return(model)
}

calc_prior_predictive_probability <- function(model){
  
  # calculate the prior predictive probability, or p(D|M)
  
  # call signature: model = calc_prior_predictive_probability(model)
  
  library(Rmpfr)
  precision = 256
  
  bic = mpfr(model$BIC, precision)
  
  model$prior_pred_prob <- exp(-0.5*bic)
  
  return(model)
  
}


calc_BF01 <- function(model_1, model_2, subject){
  
  # compute the Bayes Factor for a model comparison
  
  # evidence is calculated for first model specified in args
  
  # call signature: BF = calc_BF01(model_1, model_2)
  
  comparison = paste(model_1$Call[2], ":", model_2$Call[2])
  
  
  delta_BIC = model_2$BIC - model_1$BIC
  delta_BIC_it = model_2$BIC_it - model_1$BIC_it
  
  BF = exp(delta_BIC / 2)
  BF_it = exp(delta_BIC_it / 2)
  
  
  alt_delta_BIC = model_1$BIC - model_2$BIC
  alt_delta_BIC_it = model_1$BIC_it - model_2$BIC_it
  
  BF_alt = exp(alt_delta_BIC / 2)
  BF_alt_it = exp(alt_delta_BIC_it / 2)
  
  
  BF = list('BF'=BF, 'BF_it'=BF_it,
            'comparison'=comparison, 'BF_alt'=BF_alt,
            'BF_alt_it'=BF_alt_it, 'subject'=subject)
  
  return(BF)
}



calc_posterior_prob_h <- function(BF){
  
  # compute the posterior probability of some hypothesis given the data 
  # assuming equal prior prob. of each model
  
  # call signature: BF = calc_posterior_prob_h(BF)
  
  
  # p(H_i | data)
  BF$posterior_prob = BF$BF / (1 + BF$BF) 
  BF$posterior_prob_it = BF$BF_it / (1 + BF$BF_it)
  
  # p(H_ii | data)
  BF$posterior_prob_alt = 1 - BF$posterior_prob
  BF$posterior_prob_alt_it = 1 - BF$posterior_prob_it
  
  return(BF)
}

calc_posterior_prob_mult_h <- function(source_model, all_models){
  
  
  # assume that the prior pred probs are alrady calculated 
  
  prior_pred_prob_list = vector('list', length=length(all_models))
  
  for (model_n in seq(length(all_models))){
    
    prior_pred_prob_list[model_n] = all_models[[model_n]]$prior_pred_prob
  } 
  
  all_prior_pred_prob_sum = sum(unlist(prior_pred_prob_list))
  
  # calculate ratio
  posterior_prob_source_model = source_model$prior_pred_prob / all_prior_pred_prob_sum
  posterior_prob_source_model = posterior_prob_source_model[[1]]
  
  return(posterior_prob_source_model)
}



calc_empirical_ci <- function(BF, confidence=.95){
  
  # call signature: BF = calc_empirical_ci(BF)
  
  BF$BF_mu <- mean(BF$BF_it)
  BF$posterior_prob_mu <- mean(BF$posterior_prob_it)
  
  conf_tail = (1 - confidence) / 2
  conf_percentiles = c(conf_tail, confidence + conf_tail)
  
  BF$BF_ci = quantile(BF$BF_it, conf_percentiles)
  BF$posterior_prob_ci = quantile(BF$posterior_prob_it, conf_percentiles)
  
  return(BF)
}

# calc_empirical_ci <- function(BF, confidence=.95){
#   
#   # call signature: BF = calc_empirical_ci(BF)
#   
#   BF$BF_mu <- mean(BF$BF_it)
#   BF$posterior_prob_mu <- mean(BF$posterior_prob_it)
#   
#   conf_tail = (1 - confidence) / 2
#   conf_percentiles = c(conf_tail, confidence + conf_tail)
#   
#   BF$BF_ci = quantile(BF$BF_it, conf_percentiles)
#   BF$posterior_prob_ci = quantile(BF$posterior_prob_it, conf_percentiles)
#   
#   return(BF)
# }


calc_empirical_ci_BIC <- function(model, confidence=.95){
  
  # call signature: model = calc_empirical_ci_BIC(model)
  
  
  model$BIC_mu <- mean(model$BIC_it)
  
  conf_tail = (1 - confidence) / 2
  conf_percentiles = c(conf_tail, confidence + conf_tail)
  
  model$BIC_ci = quantile(model$BIC_it, conf_percentiles)
  
  return(model)
  
}




for (subject in seq(nrow(pc_models))){
  for (model in colnames(pc_models)[-1]){
      
      pc_models[[model]][[subject]] = calc_BIC(pc_models[[model]][[subject]])    
      pc_models[[model]][[subject]] = calc_prior_predictive_probability(pc_models[[model]][[subject]])
        
  }}

# calculating the BF for all alt models against the null 

null_BFs <- vector("list", (ncol(pc_models)-2)*4)

i = 0

for (subject in seq(nrow(pc_models))) {
  for (model in colnames(pc_models)[c(-1, -2)]) {
    # all except the subj_id and absolute null model
    
    i = i + 1
    
    null_BFs[[i]] = calc_BF01(pc_models[[model]][[subject]],
                              pc_models[[2]][[subject]],
                              pc_models$subj_id[subject]) # alt v null model
  }
}


null_BF_objs <- vector("list", (ncol(pc_models)-2)*4)

i = 0

for (subject in seq(nrow(pc_models))) {
  for (model in colnames(pc_models)[c(-1, -2)]) {
    # all except the subj_id and absolute null model
    
    i = i + 1
    
    null_BF_objs[[i]] = calc_BF01(pc_models[[model]][[subject]],
                              pc_models[[2]][[subject]],
                              pc_models$subj_id[subject]) # alt v null model
    print(paste(null_BF_objs[[i]]$subject, null_BF_objs[[i]]$comparison))
  }
}


BFs = vector("list", (ncol(pc_models)-2)*4)
subj_id = vector("list", (ncol(pc_models)-2)*4)
BF_alt = vector("list", (ncol(pc_models)-2)*4)
comparison = vector("list", (ncol(pc_models)-2)*4)

i = 0

  for (obj in null_BF_objs) {

    i = i + 1
  
    subj_id[[i]] = null_BFs[[i]]$subject
    BFs[[i]] = null_BFs[[i]]$BF
    BF_alt[[i]] = null_BFs[[i]]$BF_alt
    comparison[[i]] = null_BFs[[i]]$comparison
    
  }

pc_alt_v_null_df <- data.frame(BF=matrix(unlist(BFs), nrow=20, byrow=T),
                 subj_id=matrix(unlist(subj_id), nrow=20, byrow=T), BF_alt=matrix(unlist(BF_alt), nrow=20, byrow=T), comparison=matrix(unlist(comparison), nrow=20, byrow=T), 
stringsAsFactors=FALSE)

pc_alt_v_null_df$comparison <- as.factor(pc_alt_v_null_df$comparison)

pc_alt_v_null_df$comparison_reordered <- factor(pc_alt_v_null_df$comparison, levels=c("theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted : theta_radians_z ~ epoch_trial",
                                                                                      "theta_radians_z ~ epoch_trial + projection_0_shifted + epoch_trial * projection_0_shifted : theta_radians_z ~ epoch_trial",
                                                                                      "theta_radians_z ~ epoch_trial + projection_1_shifted + epoch_trial * projection_1_shifted : theta_radians_z ~ epoch_trial",
                                                                                      "theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + epoch_trial * projection_0_shifted + epoch_trial * projection_1_shifted : theta_radians_z ~ epoch_trial",
                                                                                      "theta_radians_z ~ epoch_trial + projection_0_shifted + projection_1_shifted + projection_0_shifted * projection_1_shifted * epoch_trial : theta_radians_z ~ epoch_trial"))

pc_alt_v_null_df$subj_id <- as.factor(pc_alt_v_null_df$subj_id)
pc_alt_v_null_df$subj_id_anon <- recode_factor(pc_alt_v_null_df$subj_id, `786` = 1, `787` = 2, `788` = 3, '789' = 4)





# BICs 

BICs = vector("list", (ncol(pc_models)-1)*4)
models = vector("list", (ncol(pc_models)-1)*4)
subj_id = vector("list", (ncol(pc_models)-1)*4)

i = 0 

for (subject in seq(nrow(pc_models))){
  for (model in colnames(pc_models)[-1]){
    
    i = i + 1
    
    BICs[[i]] = pc_models[[model]][[subject]]$BIC
    models[[i]] = model
    subj_id[[i]] = pc_models$subj_id[[subject]]
    
  }}


bic_df <- data.frame(bic=matrix(unlist(BICs), nrow=24, byrow=T),
                               subj_id=as.factor(matrix(unlist(subj_id), nrow=24, byrow=T)), 
                     model=matrix(unlist(models), nrow=24, byrow=T), 
                               stringsAsFactors=FALSE)

bic_df <- bic_df %>% group_by(subj_id) %>% mutate(null_adj_bic = bic - bic[model=='null_m'])


bic_df$subj_id_anon <- recode_factor(bic_df$subj_id, `786` = 1, `787` = 2, `788` = 3, '789' = 4)


bic_df$model <- as.factor(bic_df$model)

bic_df$model_reordered <- factor(bic_df$model, levels=c("null_m", 
                                                        "null_pc_interaction_m",
                                                        "pc0_time_interaction_m", 
                                                        "pc1_time_interaction_m", 
                                                        "pc0time_pc1time_interaction_m",
                                                        "pc01time_interaction_m"))


                                                                          
# bic_df to csv
write.csv(bic_df, file=paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_bic_df.csv'),row.names=FALSE)

# pc_alt_v_null_df to csv 
write.csv(pc_alt_v_null_df, file=paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_BF_df.csv'), row.names=FALSE)


calc_posterior_prob_mult_h <- function(source_model_prior_pred_prob, all_models_prior_pred_prob){
  
  
  # assume that the prior pred probs are alrady calculated 

  all_prior_pred_prob_sum = Reduce("+", all_models_prior_pred_prob)

  # # calculate ratio
  posterior_prob_source_model = asNumeric(source_model_prior_pred_prob[[1]] / all_prior_pred_prob_sum)

  return(posterior_prob_source_model)
}


prior_pred_prob = vector("list", (ncol(pc_models)-1)*4)
models = vector("list", (ncol(pc_models)-1)*4)
subj_id = vector("list", (ncol(pc_models)-1)*4)

i = 0 

for (subject in seq(nrow(pc_models))){
  for (model in colnames(pc_models)[-1]){
    
    i = i + 1
    
    print(paste0(subject))
    
    prior_pred_prob[[i]] <- pc_models[[model]][[subject]]$prior_pred_prob
    
    print(prior_pred_prob)
    
    models[[i]] = model
    subj_id[[i]] = pc_models$subj_id[[subject]]
    
  }
  }
  
prior_pred_prob_df <- data.frame(prior_pred_prob=matrix(prior_pred_prob, nrow=24, byrow=T),
                               subj_id=matrix(unlist(subj_id), nrow=24, byrow=T), model=matrix(unlist(models), nrow=24, byrow=T), 
                               stringsAsFactors=FALSE)

prior_pred_prob_df$subj_id_anon <- recode_factor(prior_pred_prob_df$subj_id, `786` = 1, `787` = 2, `788` = 3, '789' = 4)


subj_id = vector("list", (nrow(pc_models)))
posterior_prob_source_model = vector("list", (nrow(pc_models)))

i = 0 


for (subject in unique(prior_pred_prob_df$subj_id_anon)){
  
  i = i + 1
  
  sub_data = subset(prior_pred_prob_df, subj_id_anon == subject)
  source_model_prior_pred_prob = sub_data[which(sub_data$model == 'null_m'),]$prior_pred_prob
  all_models_prior_pred_prob = sub_data$prior_pred_prob
  # 
  # print(source_model_prior_pred_prob)
  # print(all_models_prior_pred_prob)
  
  posterior_prob_source_model[[i]] = calc_posterior_prob_mult_h(source_model_prior_pred_prob, all_models_prior_pred_prob)
  subj_id[[i]] = subject
}

posterior_prob_null_df <- data.frame(posterior_prob_null=matrix(unlist(posterior_prob_source_model), nrow=4, byrow=T),
                                 subj_id=matrix(unlist(subj_id), nrow=4, byrow=T), 
                                 stringsAsFactors=FALSE)


# posterior_prob_null_df to csv
write.csv(posterior_prob_null_df, file=paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_postprob_null_m.csv'),row.names=FALSE)

# calc posterior prob of each model P(M|D)
n_models = ncol(pc_models)-1
subj_id = vector("list", (nrow(pc_models))*n_models)
models = vector("list", (nrow(pc_models))*n_models)
posterior_prob_all_models = vector("list", (nrow(pc_models))*n_models)

i = 0 


for (subject in unique(prior_pred_prob_df$subj_id_anon)){
  
  
  sub_data = subset(prior_pred_prob_df, subj_id_anon == subject)
  
  for (model in colnames(pc_models)[-1]){
    
    i = i + 1
    
    
    print(paste('calculating p(m|d) for ', model, 'subject ', subject))
    
    source_model_prior_pred_prob = sub_data[which(sub_data$model == model),]$prior_pred_prob 
    all_models_prior_pred_prob = sub_data$prior_pred_prob
    
    posterior_prob_all_models[[i]] = calc_posterior_prob_mult_h(source_model_prior_pred_prob, all_models_prior_pred_prob)
    subj_id[[i]] = subject
    models[[i]] = model
    

  }
 
}

all_models_posterior_prob_df <- data.frame(posterior_prob=matrix(unlist(posterior_prob_all_models), nrow=24, byrow=T),
                                     subj_id=factor(matrix(unlist(subj_id), nrow=24, byrow=T)), model=factor(matrix(unlist(models), nrow=24, byrow=T)), 
                                     stringsAsFactors=FALSE)

fig_dir <- paste0(home, 'Dropbox/loki_0.5/figures/model_comparisons/')



all_models_posterior_prob_df$model_reordered <- factor(all_models_posterior_prob_df$model, 
                                                       levels=c("null_m", 
                                                        "null_pc_interaction_m",
                                                        "pc0_time_interaction_m", 
                                                        "pc1_time_interaction_m", 
                                                        "pc0time_pc1time_interaction_m",
                                                        "pc01time_interaction_m"))



# all_models_posterior_prob_df to csv
write.csv(all_models_posterior_prob_df, file=paste0(home, 'Dropbox/loki_0.5/analysis/aggregated_data/pc_postprob_each_m.csv'),row.names=FALSE)

