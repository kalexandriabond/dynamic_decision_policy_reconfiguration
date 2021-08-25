import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import hddm
import numpy as np
from sys import platform
import os
# print(hddm.__version__)

# %matplotlib inline
# plt.rcParams['figure.figsize'] = 30, 10
sns.set_context("notebook", font_scale=2)


if platform == 'linux2':
    home = '/lab_data/coaxlab/Projects/dynamic_decision_policy_reconfiguration/'
elif platform == 'darwin':
    home = '/Users/i_67981492/dynamic_decision_policy_reconfiguration/'

print(platform, home)

write_dir = os.path.join(home, 'analysis_revision/loki0_stim_avtz_dc_models'); print(write_dir)

os.path.isdir(write_dir)

all_obs_data = hddm.load_csv(os.path.join(home, 'aggregated_data/loki_0/av_est_fix.csv'))
all_obs_data_pared = all_obs_data.loc[(all_obs_data.shifted_epoch_trial <= 3) &
                                      (all_obs_data.shifted_epoch_trial >=-1)]

all_obs_data_pared.groupby('condition').p_optimal.unique()

all_obs_data_pared.groupby('condition').lambda_val.unique()


all_obs_data_pared.columns = all_obs_data_pared.columns.str.strip()
all_obs_data_pared.head()


all_obs_data_pared.loc[all_obs_data_pared.cue_choice == 111, 'id_choice'] = 0
all_obs_data_pared.loc[all_obs_data_pared.cue_choice == 112, 'id_choice'] = 1

all_obs_data_pared.loc[all_obs_data_pared.high_p_cue == 111, 'p_id_solution'] = 0
all_obs_data_pared.loc[all_obs_data_pared.high_p_cue == 112, 'p_id_solution'] = 1


all_obs_data_pared["shifted_epoch_trial"] = all_obs_data_pared["shifted_epoch_trial"].astype('category')

all_obs_data_pared_stim = all_obs_data_pared.rename(index=str, columns={"id_choice": "response",
                                                       "p_id_solution": "stim",
                                                        "subj_id": "subj_idx"})
all_obs_data_pared_stim = all_obs_data_pared_stim[['response', 'stim', 'rt', 'experiment',
                                                   'condition',
                             'subj_idx', 'shifted_epoch_trial']].dropna()

included_stim_params_input = ['a', 'v', 't', 'z']
included_stim_params = ['a', 'v', 't', 'z', 'dc']

split_param = 'v'

n_samples, n_burned_samples, n_thin = 20000, 5000, 5
p_outlier = 0.05
#n_samples, n_burned_samples, n_thin = 50, 5, 2 # min. samples to test the workflow


# estimate models

null_stim_m = hddm.HDDMStimCoding(all_obs_data_pared_stim, stim_col='stim',
              include=included_stim_params_input,
                        p_outlier=p_outlier,
                        trace_subjs=True,
                        informative=True,
                                  group_only_nodes = included_stim_params,
                                 drift_criterion=True,
                                 split_param=split_param)
null_stim_m.find_starting_values()
null_stim_m.sample(n_samples, burn=n_burned_samples, thin=n_thin)


models = []

for param in included_stim_params:

    depends_on = {param: 'shifted_epoch_trial'}
    print(depends_on)

    m = hddm.HDDMStimCoding(all_obs_data_pared_stim, stim_col='stim',
                            depends_on=depends_on,
                  include=included_stim_params_input,
                            p_outlier=p_outlier,
                            trace_subjs=True,
                            informative=True,
                            group_only_nodes = included_stim_params,
                           drift_criterion=True,
                           split_param=split_param)
    m.find_starting_values()
    m.sample(n_samples, burn=n_burned_samples, thin=n_thin)

    models.append(m)

# plot posteriors to assess convergence

for m in models:

    plt.figure()
    m.plot_posteriors(save=True)
    plt.title(m.depends_on)

# save dics

model_dic_list = []

null_dic_dict = {'model': 'null'}
null_dic_dict.update(null_stim_m.dic_info)

model_dic_list.append(null_dic_dict)

for m in models: # alternatives

    depends_on_var = m.depends_on.values()[0][0]
    ddm_param = m.depends_on.keys()[0]

    temp_dict = {'model': ddm_param + '_' + depends_on_var}
    temp_dict.update(m.dic_info)

    temp_dict['null_adj_DIC'] = null_dic_dict['DIC'] - temp_dict['DIC']

    model_dic_list.append(temp_dict)


dic_df = pd.DataFrame(model_dic_list)
dic_df.to_csv(os.path.join(home, write_dir, 'avtz_dc_dic_df.csv'), index=False)


# save traces

group_trace_dfs = []

for m in models:

    group_traces = m.get_group_traces()

   # identify the model
    depends_on_var = m.depends_on.values()[0][0]
    ddm_param = m.depends_on.keys()[0]

    # group traces
    group_trace_cols = [col for col in group_traces.columns if ddm_param in col]
    group_only_traces = group_traces[group_trace_cols]

    group_only_traces['mcmc_iteration'] = np.arange(len(group_only_traces))
    #group_only_traces = group_traces.drop(columns= ddm_param + '_std')

    group_only_traces_melted = pd.melt(group_only_traces, var_name='shifted_epoch_trial', value_name=ddm_param, id_vars='mcmc_iteration')
    group_only_traces_melted['shifted_epoch_trial'] = group_only_traces_melted.shifted_epoch_trial.str.split('(', expand=True)[1].str.split(')', expand=True)[0]


    group_only_traces_melted['model'] = (ddm_param + '_' + depends_on_var)
    group_only_traces_melted['ddm_param'] = ddm_param
    group_only_traces_melted.head()

    group_only_traces_melted.to_csv(os.path.join(home, write_dir,
                                                 ddm_param + '_' + depends_on_var
                                                 + '_traces.csv'), index=False)

    group_trace_dfs.append(group_only_traces_melted)


# plot CP-evoked decision param dists

for trace_df in group_trace_dfs:

    ddm_param = trace_df.ddm_param.unique()[0]

    sns.set_context("notebook", font_scale=1.7)
    sns.set_style('white')
    sns.despine()

    g = sns.FacetGrid(trace_df,
                      palette="Set1", size=5, aspect=3)

    g = (g.map(sns.violinplot, "shifted_epoch_trial", ddm_param, split=True))
    plt.xlabel('Trials after change point')

    g.savefig(os.path.join(home, write_dir, ddm_param + '_cp_evoked_distributions.png'))
