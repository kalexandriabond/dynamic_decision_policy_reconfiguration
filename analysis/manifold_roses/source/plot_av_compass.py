import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns
import time
import os

from matplotlib import rc
# from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import circmean


from jupyterthemes import jtplot
jtplot.style(context='talk', fscale=1.1, spines=False, gridlines='--', )

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('xtick', labelsize=10)
rc('figure', figsize=(10, 5))

def rose_plot(ax, angles, bins=16, density=None, offset=0, lab_unit="degrees",
              start_zero=False, **param_dict):
    """
    Plot polar histogram of angles on ax. ax must have been created using
    subplot_kw=dict(projection='polar'). Angles are expected in radians.
    """

    clean_angles = angles[~np.isnan(angles)]
    # Wrap angles to [-pi, pi)
    angles = (clean_angles + np.pi) % (2*np.pi) - np.pi

    # Set bins symmetrically around zero
    if start_zero:
        # To have a bin edge at zero use an even number of bins
        if bins % 2:
            bins += 1
        bins = np.linspace(-np.pi, np.pi, num=bins+1)

    # Bin data and record counts
    count, bin = np.histogram(angles, bins=bins)

    # Compute width of each bin
    widths = np.diff(bin)

    # By default plot density (frequency potentially misleading)
    # Radius of histogram is proportional to the density of angular observations
    if density is None or density is True:
        # Area to assign each bin
        area = count / angles.size
        # Calculate corresponding bin radius
        radii = (area / np.pi)**.5
    else:
        radii = count

    # Plot data on ax
    ax.bar(bin[:-1], radii, zorder=1, align='edge', width=widths,
           edgecolor='gray', fill=True, linewidth=1, alpha=0.2, color='gray')

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels, they are mostly obstructive and not informative
    ax.set_yticks([])

    if lab_unit == "radians":
        label = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',
                  r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$']
        ax.set_xticklabels(label)


    return radii


def plot_roses(av_manifold_df_unit, n_rows=2, n_cols=5, n_conditions=9,
savefig=None, conditional=False, all_subs_all_conds=True, pub_quality_savefig=None, home = os.path.expanduser('~'), id_str=None, sample_sub=None):

    fig_path = os.path.join(home, 'Dropbox/loki_0_0.5/loki_0.5_figures/rose/')

    fig, main_ax = plt.subplots(n_rows, n_cols, subplot_kw=dict(polar=True))

    ax_list = main_ax.flatten()
    kw = dict(arrowstyle="->", color='gray', lw=0.5, alpha=0.3)
    mean_kw = dict(arrowstyle="->", color='black', lw=1.5, alpha=0.7)
    shifted_epoch_trials = np.sort(av_manifold_df_unit.shifted_epoch_trial.unique())[1:] # don't need first, just np.nan

    for shifted_trial, ax in zip(shifted_epoch_trials, ax_list):
        data = av_manifold_df_unit.loc[av_manifold_df_unit.shifted_epoch_trial == shifted_trial].reset_index(drop=True).copy()

        radii = rose_plot(ax, data.theta_radians_z)
        plotting_radius = radii.max()

        # radii = counts, not vec length
        # way to normalize reference frame? to plot both radii & mean? 

        if shifted_trial == 0:
            ax.set_title(r'$\vec{{_{{-\textrm{}:\textrm{}}}}}$'.format(abs(int(shifted_trial-1)), int(shifted_trial)), fontsize=15, y=1.3)
        else:
            ax.set_title(r'$\vec{{_{{\textrm{}:\textrm{}}}}}$'.format(int(shifted_trial)-1, int(shifted_trial)), fontsize=15, y=1.3)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.tick_params(axis='y', which='major', labelsize=8)
        plt.tick_params(axis='x', which='minor', labelsize=15)

        print('r_max: ', radii.max()) # these are counts for angle, not vector length
        if radii.max() is not np.nan:
            ax.set_ylim(0, radii.max())

        mean_theta = circmean(data.theta_radians_z, nan_policy='omit') # assumes radians as input
        print('mean theta ', mean_theta, 'shifted trial ', data.shifted_epoch_trial.unique(), 'subject ', data.subj_id.unique(), 'condition ', data.condition.unique())
        ax.annotate("", xy=(mean_theta, plotting_radius), xytext=(0,0), arrowprops=mean_kw)

    if (len(ax_list)%2) > 1:
        fig.delaxes(ax_list[-1]) # if odd number of plots delete last subplot

    if conditional:
        if id_str:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + str(int(av_manifold_df_unit.condition.unique()[0])) +'_av_manifold_polar_rose' + id_str + '.png')
        else:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + str(int(av_manifold_df_unit.condition.unique()[0])) +'_av_manifold_polar_rose.png')

        fig.suptitle('subject ' + str(int(av_manifold_df_unit.subj_id.unique()[0])) + ': ' +
        '$\lambda =$ ' + str(int(av_manifold_df_unit.lambda_val.unique()[0])) + ' p = ' + str(av_manifold_df_unit.p_optimal.unique()[0]), fontsize=20)
    elif sample_sub:
        if id_str:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + 'all_conditions_av_manifold_polar_rose' + id_str + '.png')
    elif all_subs_all_conds:
        if id_str:
            fig_name = ('all_subs_all_conditions_av_manifold_polar_rose' + id_str + '.png')
        else:
            fig_name = ('all_subs_all_conditions_av_manifold_polar_rose.png')

    if savefig:
        # plt.get_current_fig_manager().window.state('zoomed')
        plt.get_current_fig_manager().window.showMaximized()
        plt.tight_layout(pad=3.0) # control spacing between subplots
        plt.savefig(os.path.join(fig_path, fig_name), dpi=1200)
        plt.show()
        plt.pause(0.1)
        plt.close()



    return fig, ax


av_manifold_df = pd.read_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv')

av_manifold_subset_df = av_manifold_df.loc[(av_manifold_df.shifted_epoch_trial >= -1) & (av_manifold_df.shifted_epoch_trial < 4)]


# sanity check on random sub.

random_sample_sub = av_manifold_subset_df.subj_id.sample(n=1, random_state=11620).values[0]

random_sample_sub_data = av_manifold_subset_df.loc[(av_manifold_subset_df.subj_id == random_sample_sub)].reset_index(drop=True).copy()


n_conditions = 9

assert random_sample_sub_data.condition.nunique() == n_conditions, 'check data str'
assert random_sample_sub_data.subj_id.nunique() == 1, 'check data str'

_,_ = plot_roses(av_manifold_subset_df, savefig=True, n_rows=1, n_cols=4, all_subs_all_conds=True)
