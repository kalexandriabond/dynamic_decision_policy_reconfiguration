
library(latex2exp)
library(tidyverse)
library(circglmbayes)

library(reshape2)

load("~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/circ_glm/angular_decay.RData")

trace_df <- as.data.frame.array(angular_decay_m$all_chains)

head(trace_df)

fig_dir = '~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/figures/circ_glm_figures/'

deg_trace_df = trace_df[-2] * 180/pi # exclude kappa param from conversion

dt_chain_df <- melt(trace_df[,3:10], variable.name = 'dt_chain', value.name = 'step_rad')

mu_chain_df <- melt(trace_df[,11:19], variable.name = 'mu_chain', value.name = 'step_rad')

write.csv(row.names = FALSE, x = deg_trace_df, file = '~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/circ_glm/deg_trace_df.csv')
write.csv(row.names = FALSE, x = dt_chain_df, file = '~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/circ_glm/dt_chain_df.csv')
write.csv(row.names = FALSE, x = mu_chain_df, file = '~/Dropbox/loki_0/simple_rt_experiment_probabilityC/analysis/circ_glm/mu_chain_df.csv')


ggplot(data=trace_df, aes(x=b0_chain)) + geom_histogram() +
  xlab(TeX('$\\beta_0$'))
ggplot(data=trace_df, aes(x=kp_chain)) + geom_histogram() +
  xlab(TeX('$\\kappa$'))


mu_chain_df$step_deg = (mu_chain_df$step_rad * 180/pi)
dt_chain_df$step_deg = (dt_chain_df$step_rad * 180/pi)


ggplot(data=mu_chain_df, aes(x=mu_chain, y=step_deg)) + geom_violin(fill='gray') +
  ylab(TeX('Estimated $\\theta$ $(\\degree)$')) + xlab('Vector')  +
  scale_x_discrete(labels = c(parse(text=TeX('$\\Delta_{-1:0}$')), 
                              parse(text=TeX('$\\Delta_{0:1}$')),
                              parse(text=TeX('$\\Delta_{1:2}$')), 
                              parse(text=TeX('$\\Delta_{2:3}$')),
                              parse(text=TeX('$\\Delta_{3:4}$')),
                              parse(text=TeX('$\\Delta_{4:5}$')),
                              parse(text=TeX('$\\Delta_{5:6}$')),
                              parse(text=TeX('$\\Delta_{6:7}$')),
                              parse(text=TeX('$\\Delta_{7:8}$')))) + theme_minimal(base_size = 16,)
ggsave(paste0(fig_dir, 'estimated_vector_theta_deg_decay.png'))

ggplot(data=mu_chain_df, aes(x=mu_chain, y=step_rad)) + geom_violin(fill='gray') +
  ylab(TeX('Estimated $\\theta$ (rad.)')) + xlab('Vector') +
  scale_x_discrete(labels = c(parse(text=TeX('$\\Delta_{-1:0}$')), 
                              parse(text=TeX('$\\Delta_{0:1}$')),
                              parse(text=TeX('$\\Delta_{1:2}$')), 
                              parse(text=TeX('$\\Delta_{2:3}$')),
                              parse(text=TeX('$\\Delta_{3:4}$')),
                              parse(text=TeX('$\\Delta_{4:5}$')),
                              parse(text=TeX('$\\Delta_{5:6}$')),
                              parse(text=TeX('$\\Delta_{6:7}$')),
                              parse(text=TeX('$\\Delta_{7:8}$')))) + 
  theme_minimal(base_size = 16)
ggsave(paste0(fig_dir, 'estimated_vector_theta_rad_decay.png'))



ggplot(data=dt_chain_df, aes(x=dt_chain, y=step_deg)) + geom_violin(fill='gray') +
  ylab(TeX('Estimated deviation from initial vector ($\\theta\\Delta_{n} - \\theta\\Delta_{0:1}  (\\degree))$')) + xlab('Vector')  +
  
  geom_hline(yintercept = 0, 
             alpha=1, color='black', size=.8,  linetype="dashed") +

  scale_x_discrete(labels = c(
                              parse(text=TeX('$\\Delta_{0:1}$')),
                              parse(text=TeX('$\\Delta_{1:2}$')), 
                              parse(text=TeX('$\\Delta_{2:3}$')),
                              parse(text=TeX('$\\Delta_{3:4}$')),
                              parse(text=TeX('$\\Delta_{4:5}$')),
                              parse(text=TeX('$\\Delta_{5:6}$')),
                              parse(text=TeX('$\\Delta_{6:7}$')),
                              parse(text=TeX('$\\Delta_{7:8}$')))) + theme_minimal(base_size = 12,)
  
  ggsave(paste0(fig_dir, 'estimated_dt_vector_theta_deg_decay.png'))


 
BF_equal <- BF.circGLM(angular_decay_m)$PMP_Mean_Eq

sequential_p_equal = as.data.frame(c(BF_equal["[Reference, epoch_trial1]", 1],
BF_equal["[epoch_trial1, epoch_trial2]",1],
BF_equal["[epoch_trial2, epoch_trial3]",1],
BF_equal["[epoch_trial3, epoch_trial4]",1],
BF_equal["[epoch_trial4, epoch_trial5]",1],
BF_equal["[epoch_trial5, epoch_trial6]",1],
BF_equal["[epoch_trial6, epoch_trial7]",1],
BF_equal["[epoch_trial7, epoch_trial8]",1]))

colnames(sequential_p_equal) <- c('p_equal')
sequential_p_equal$dt <- c('01', '12', '23', '34', '45','56', '67', '78')


ggplot(data=sequential_p_equal, aes(x=dt, y=p_equal, group=1)) + geom_line() + geom_point() +
  ylab(TeX('posterior p($\\theta\\Delta_{n} == \\theta\\Delta_{n+1}$)')) + 
  
  scale_x_discrete(labels = c(parse(text=TeX('$\\theta\\Delta_{01}$')), 
                              parse(text=TeX('$\\theta\\Delta_{12}$')),
                              parse(text=TeX('$\\theta\\Delta_{23}$')), 
                              parse(text=TeX('$\\theta\\Delta_{34}$')),
                              parse(text=TeX('$\\theta\\Delta_{45}$')),
                              parse(text=TeX('$\\theta\\Delta_{56}$')),
                              parse(text=TeX('$\\theta\\Delta_{67}$')),
                              parse(text=TeX('$\\theta\\Delta_{78}$')))) + theme_minimal(base_size = 12,) + 
  xlab(TeX('$\\theta$ distribution comparison$')) 
ggsave(paste0(fig_dir, 'p_equal_sequential_theta_distributions.png'))

  
  


