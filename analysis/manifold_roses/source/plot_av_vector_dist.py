
import numpy as np
import pandas as pd
import seaborn as sns
import time
import os

import matplotlib.pyplot as plt

from matplotlib import rc
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


from jupyterthemes import jtplot
jtplot.style(context='talk', fscale=2, spines=False, gridlines='--', )


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('figure', figsize=(10,5))

import matplotlib.pyplot as plt


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



def plot_compass(av_manifold_df_unit, n_rows=2, n_cols=2,
n_conditions=4, plot_mean_components=None, savefig=None,
all_conditions=False, r_constant=None, id_str=None, pub_quality_savefig=None, all_data=False,
 home=os.path.expanduser('~')):

    import matplotlib.pyplot as plt

    fig_path = os.path.join(home, 'Dropbox/loki_0_0.5/loki_0.5_figures/rose/')

    fig, main_ax = plt.subplots(n_rows, n_cols, subplot_kw=dict(polar=True),
     gridspec_kw={'hspace': 0.6, 'wspace': .5, 'top': 0.88,
    'left': 0.11, 'right': .9, 'bottom': .11})


    ax_list = main_ax.flatten()
    print('ax_list len ', len(ax_list))
    kw = dict(arrowstyle="->", color='gray', lw=2, alpha=0.3)
    mean_kw = dict(arrowstyle="->", color='black', lw=2.5)

    shifted_epoch_trials = np.sort(av_manifold_df_unit.shifted_epoch_trial.unique())[1:] # don't need first, just np.nan



    for shifted_trial, ax in zip(shifted_epoch_trials, ax_list):
        data = av_manifold_df_unit.loc[av_manifold_df_unit.shifted_epoch_trial == shifted_trial].reset_index(drop=True).copy()
        r_max = av_manifold_df_unit.r_z.max()

        if r_constant:
            constant_r = r_max + np.zeros_like(data.r_z)
            [ax.annotate("", xy=(theta, radius), xytext=(0,0), arrowprops=kw) for radius, theta in zip(constant_r, data.theta_radians_z)]

            mean_rad = r_max

        else:
            [ax.annotate("", xy=(theta, radius), xytext=(0,0), arrowprops=kw) for radius, theta in zip(data.r_z, data.theta_radians_z)]
            mean_rad = np.mean(data.r_z)


        ax.yaxis.set_major_formatter(FormatStrFormatter('%1d'))
        ax.tick_params(axis='y', which='major', labelsize=20)
        ax.tick_params(axis='x', which='major', labelsize=20)

        print('r_max: ', r_max)
        if r_max is not np.nan:
            ax.set_ylim(0, r_max)


        from scipy.stats import circmean

        mean_theta = circmean(data.theta_radians_z, nan_policy='omit') # assumes radians as input
        print('mean theta ', mean_theta, 'mean rad ', mean_rad, 'shifted trial ', data.shifted_epoch_trial.unique(), 'subject ', data.subj_id.unique(), 'condition ', data.condition.unique())

        if plot_mean_components:
            rads = np.arange(0, (2*np.pi), 0.01)
            mean_theta_arr = np.ones_like(rads)*mean_theta
            plt.polar(mean_theta, r_max,  '.', color='black', clip_on=False, markersize=5)
            plt.polar(rads, mean_theta_arr, '-', color='black')

        ax.annotate("", xy=(mean_theta, mean_rad), xytext=(0,0), arrowprops=mean_kw)

        if shifted_trial == 0:
            ax.set_title(r'$\vec{{_{{-\textrm{}:\textrm{}}}}}$'.format(abs(int(shifted_trial-1)), int(shifted_trial)), fontsize=25, y=1.20)
        else:
            ax.set_title(r'$\vec{{_{{\textrm{}:\textrm{}}}}}$'.format(int(shifted_trial)-1, int(shifted_trial)), fontsize=25, y=1.20)

    if (len(ax_list)%2) > 1:
        fig.delaxes(ax_list[-1]) # odd number of plots so delete last subplot

    if all_conditions:
        if id_str:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + 'all_conditions_av_manifold_polar' + id_str + '.png')
        else:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + 'all_conditions_av_manifold_polar.png')

        fig.suptitle('subject ' + str(int(av_manifold_df_unit.subj_id.unique()[0])) +  ': all data', fontsize=20)
    if all_data:
        fig_name = ('all_conditions_all_subs_ all_subs_all_conditions_av_manifold_polar_vector_dist.png')


    else:
        if id_str:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + str(av_manifold_df_unit.condition.unique()[0]) +'_av_manifold_polar' + id_str + '.png')
        else:
            fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + str(av_manifold_df_unit.condition.unique()[0]) +'_av_manifold_polar.png')

        fig.suptitle('subject ' + str(int(av_manifold_df_unit.subj_id.unique()[0])) + ': ' +
        '$\lambda =$ ' + str(int(av_manifold_df_unit.lambda_val.unique()[0])) + ' p = ' + str(av_manifold_df_unit.p_optimal.unique()[0]), fontsize=20)


    if savefig:

        plt.get_current_fig_manager().window.showMaximized()
        plt.savefig(os.path.join(fig_path, fig_name), dpi=1200)
        plt.show()
        plt.pause(0.1)
        plt.close()

    return fig, ax_list

home = '/home/krista/'

av_manifold_df = pd.read_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv')

random_subj_id = av_manifold_df.subj_id.sample().values[0]

av_manifold_df_unit = av_manifold_df.loc[(av_manifold_df.subj_id == random_subj_id)].reset_index(drop=True).copy()


# testing with my data

test_theta = av_manifold_df_unit.theta_radians_z.values
fig, ax = plt.subplots(1,1, subplot_kw=dict(projection='polar'))
rose_plot(ax, test_theta)
kw = dict(arrowstyle="->", color='gray', lw=2, alpha=0.3)
ax.annotate("", xy=(1, 1), xytext=(0,0), arrowprops=kw) # ylim diff? probably
fig.tight_layout(); fig.show();


av_manifold_subset_df = av_manifold_df.loc[(av_manifold_df.shifted_epoch_trial >= -1) & (av_manifold_df.shifted_epoch_trial <= 4)]

random_sample_sub = av_manifold_subset_df.subj_id.sample(n=1, random_state=11620).values[0]

random_sample_sub_data = av_manifold_subset_df.loc[(av_manifold_subset_df.subj_id == random_sample_sub)].reset_index(drop=True).copy()




jtplot.style(context='talk', fscale=1.5, spines=False, gridlines='--', )

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.ioff()

_,_ = plot_compass(av_manifold_subset_df, savefig=True, n_rows=2, n_cols=2, all_data=True)
