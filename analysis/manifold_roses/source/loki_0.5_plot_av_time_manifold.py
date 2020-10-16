import numpy as np
import pandas as pd
import seaborn as sns
import time
import os
import matplotlib.pyplot as plt
from scipy import stats
import itertools

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

from jupyterthemes import jtplot
jtplot.style(context='talk', fscale=3, spines=False, gridlines='--', )

av_manifold_df = pd.read_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv')

write_dir = os.path.join(os.path.expanduser('~'), 'Dropbox/loki_0_0.5/loki_0.5_figures/')

def plot_a_v_time(data, fig, ax,  conditional=False, savefig=None,
all_subs=False, linestyle='-', legend=True,
pooled_mean=False, id_str=None, sample_sub=False, home=os.path.expanduser('~')):


    # hack to get hue to work with lineplots and markers ...

    fig_path=write_dir

    n_plotted_trials = data.shifted_epoch_trial.nunique()

    palette_seed = sns.color_palette('Greens', n_colors=200)[80::10]

    assert len(palette_seed) >= n_plotted_trials, 'check n_colors for color palette'

    palette = itertools.cycle(palette_seed)

    sns.lineplot(data=data, x='a_est_z', y='v_est_z', hue='shifted_epoch_trial', palette=palette_seed[:n_plotted_trials], marker='o');

    x = data.a_est_z
    y = data.v_est_z

    for i in range(len(data)):
        plt.plot(x.values[i:i+2], y.values[i:i+2], color=next(palette), linestyle=linestyle, linewidth=4)

    if conditional is True:
        plt.title('subject ' + str(int(data.subj_id.unique()[0])) + ': ' +
                '$\lambda =$ ' + str(int(data.lambda_val.unique()[0])) + ' p = ' + str(data.p_optimal.unique()[0]), fontsize=20)
        fig_name = (id_str + str(int(data.subj_id.unique()[0])) + '_' + str(int(data.condition.unique()[0])) +'_a_v_time_color.png')
    if all_subs is True:
        fig_name = (id_str + 'all_conditions_all_subs_a_v_time_color.png')
    if sample_sub is True:
        fig_name = (id_str + 'sample_sub_all_conditions_a_v_time_color.png')
    if pooled_mean is True:
        fig_name = (id_str + 'mean_a_v_time_color.png')
    # else:
    #     # plt.title('subject ' + str(int(data.subj_id.unique()[0])), fontsize=20)
    #     fig_name = (str(int(data.subj_id.unique()[0])) +'_all_conditions_a_v_time_color.png')
    #

    plt.xlabel(r'Boundary height ($\hat{a}$)')
    plt.ylabel(r'Drift rate ($\hat{v}$)')

    if legend:
        legend = ax.legend()
        legend.texts[0].set_text("epoch trial")
    else:
        ax.get_legend().remove()


    if savefig:
        plt.tight_layout() # control spacing between subplots
        plt.savefig(os.path.join(fig_path, fig_name))


    return fig, ax


av_manifold_df_subset = av_manifold_df.loc[(av_manifold_df.shifted_epoch_trial >= -1) & (av_manifold_df.shifted_epoch_trial < 4)].reset_index(drop=True).copy()

mean_av_df = av_manifold_df_subset.groupby(['shifted_epoch_trial'])[['a_est_z', 'v_est_z']].mean().reset_index() # only plot t0-4

sub_mean_av_df = av_manifold_df_subset.groupby(['subj_id', 'shifted_epoch_trial'])[['a_est_z', 'v_est_z']].mean().reset_index() # only plot t0-6


# mean
fig, ax = plt.subplots()
plot_a_v_time(mean_av_df, fig, ax, conditional=False, savefig=True,
all_subs=True, pooled_mean=True, legend=False, id_str='loki_0.5_')

# sample sub
random_sub = sub_mean_av_df.subj_id.sample(random_state=92020).values[0]
sub_data = sub_mean_av_df.loc[sub_mean_av_df.subj_id == random_sub].reset_index().copy()

fig, ax = plt.subplots()
plot_a_v_time(sub_data, fig, ax, conditional=False, savefig=True, all_subs=False, sample_sub=True,
legend=False, id_str='loki_0.5_')
plt.show()

# all subs
fig, ax = plt.subplots()

for subj_id in sub_mean_av_df.subj_id.unique():

    sub_data = sub_mean_av_df.loc[sub_mean_av_df.subj_id == subj_id].reset_index().copy()

    plot_a_v_time(sub_data, fig, ax, conditional=False, savefig=True, all_subs=True, legend=False, id_str='loki_0.5_')
