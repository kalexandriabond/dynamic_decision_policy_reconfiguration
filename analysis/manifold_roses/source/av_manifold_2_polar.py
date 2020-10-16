import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import time
import os

av_manifold_df = pd.read_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv')

assert 'a_est_z' in av_manifold_df, 'z-scored a is missing. check df.'
assert 'v_est_z' in av_manifold_df, 'z-scored v is missing. check df.'


# testing the data
# assumption: theta and r for first shifted_epoch_trial of each epoch window number must be np.nan (single point, diff is np.nan)
first_trial_epoch_window_df = av_manifold_df.groupby(['subj_id', 'condition', 'epoch_window_number'])['theta_radians', 'r'].apply(lambda x: x.head(1)).reset_index()

assert (first_trial_epoch_window_df.theta_radians.isna()).sum() == len(first_trial_epoch_window_df), 'check first trial of epoch window'
assert (first_trial_epoch_window_df.r.isna()).sum() == len(first_trial_epoch_window_df), 'check first trial of epoch window'


# assumption 2: max shifted_epoch_trial in epoch_window_number should NOT be negative and should correspond to the length of the epoch
# the tail of shifted_epoch_trial and epoch_trial should be the same except for epoch_window_number 0


# shifted_epoch_trial lsat one is not always 8

epoch_trial_shifted_epoch_trial_df = av_manifold_df.groupby(['subj_id', 'condition', 'epoch_window_number'])[['epoch_trial', 'shifted_epoch_trial']].apply(lambda x: x.tail(1)).reset_index()






def cart2polar(x,y): # note that this should be implemented for each epoch after nesting within condition and sub

    x_diff, y_diff = x.diff(), y.diff()
    r = np.hypot(x_diff, y_diff)
    theta_radians = np.arctan2(y_diff, x_diff)
    theta_deg = np.rad2deg(theta_radians)
    return r, theta_radians, theta_deg


# shfited_epoch_trial is not always as long as the desired window of -1:8
def create_epoch_window_number(session_data):


    print('processing sub {} session {}'.format(session_data.subj_id.unique(), session_data.condition.unique()))
    session_data['epoch_window_number'] = session_data.epoch_number

    for epoch_number in session_data.epoch_number.unique()[1:]:  # only get the second epoch onward (otherwise no prior trial to include)
        epoch_start_idx = session_data.loc[session_data.epoch_number == epoch_number].index[0]
        epoch_window_start_idx = epoch_start_idx - 1

        print('epoch start idx ', epoch_start_idx)
        print('epoch_window_start_idx ', epoch_window_start_idx)

        print('epoch_start_idx > epoch_window_start_idx ', epoch_start_idx > epoch_window_start_idx)

        print('epoch start idx ', epoch_start_idx)
        print('epoch window start idx ', epoch_window_start_idx)

        session_data.set_value(epoch_window_start_idx, 'epoch_window_number',
        session_data.iloc[epoch_start_idx].epoch_number)

    return session_data

reg_window_dfs = []

# now implement on all data
for subj_id in av_manifold_df.subj_id.unique():
    for condition in av_manifold_df.condition.unique():
        session_data = av_manifold_df.loc[(av_manifold_df.subj_id == subj_id) &
        (av_manifold_df.condition == condition)].reset_index(drop=True).copy()

        session_data_v2 = create_epoch_window_number(session_data)

        # print(session_data_v2[['epoch_window_number', 'epoch_number', 'shifted_epoch_trial', 'epoch_trial']].head(60))


        reg_window_dfs.append(session_data_v2)


reg_window_df_v2 = pd.concat(reg_window_dfs, axis=0).reset_index(drop=True)



epoch_datumz = []

for subj_id in av_manifold_df.subj_id.unique():
    for condition in av_manifold_df.condition.unique():

        session_data = av_manifold_df.loc[(av_manifold_df.subj_id == subj_id) & (av_manifold_df.condition == condition)].copy()
        for epoch_window_number in session_data.epoch_window_number.unique(): # how is epoch_window_number being calculated? including data from multiple epochs?

            epoch_data = session_data.loc[session_data.epoch_window_number == epoch_window_number].reset_index(drop=True).copy()

            r, theta_radians, theta_deg = cart2polar(epoch_data.a_est_z, epoch_data.v_est_z)

            epoch_data['r_z'] = r
            epoch_data['theta_radians_z'] = theta_radians
            epoch_data['theta_deg_z'] = theta_deg

            epoch_datumz.append(epoch_data)

av_manifold_df = pd.concat(epoch_datumz, axis=0)
av_manifold_df.to_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv', index=False)

g = sns.FacetGrid(data=av_manifold_df, col='subj_id')
g = g.map(plt.hist, 'theta_radians_z')
plt.savefig('theta_radians_z.png')

g = sns.FacetGrid(data=av_manifold_df, col='subj_id')
g = g.map(plt.hist, 'r_z')
plt.savefig('r_z.png')

g = sns.FacetGrid(data=av_manifold_df, col='subj_id')
g = g.map(plt.hist, 'theta_deg_z')
plt.savefig('theta_deg_z.png')

# ___

def plot_compass(r_z, theta_radians_z, epoch_trial_n,  arrowprops=None, r_max=1,
plot_mean_components=False,  n_rows=1, n_cols=8, color='black'): # note that this takes radians but plots in degrees
    fig, ax = plt.subplots(n_rows, n_cols, subplot_kw=dict(polar=True))
    ax_list = ax.flatten()

    kw = dict(arrowstyle="->", color='gray', lw=2)
    mean_kw = dict(arrowstyle="->", color='black', lw=2.5)

    for ax in ax_list:
        [ax.annotate("", xy=(theta, rad), xytext=(0,0), arrowprops=kw) for rad, theta in zip(r_z, theta_radians_z)]
        ax.set_ylim(0, r_max)



    mean_rad = np.mean(r_z)
    mean_theta = np.mean(theta_radians_z)

    if plot_mean_components:

        rads = np.arange(0, (2*np.pi), 0.01)
        mean_theta_arr = np.ones_like(rads)*mean_theta

        plt.polar(mean_theta, r_max,  '.', color='black', clip_on=False, markersize=20)

        plt.polar(rads, mean_theta_arr, '-', color='black')

    ax.annotate("", xy=(mean_theta, mean_rad), xytext=(0,0), arrowprops=mean_kw)

    plt.title('epoch trial vector '+ str(epoch_trial_n-1) + ':' + str(epoch_trial_n), y=1.08)

    return fig, ax


def evaluate_magnitude_uniformity(trial_data):
    # from astropy.stats import rayleightest ?

    # plot distribution
    plt.figure()
    plt.title('distribution of magnitudes for epoch trial vector {}:{}'.format((trial_data.shifted_epoch_trial.unique()-1)[0],
    trial_data.shifted_epoch_trial.unique()[0]))
    plt.hist(trial_data.r)
    plt.xlabel('magnitude (radius)'); plt.ylabel('frequency')


    # evaluate statistical sig.
    # p_uniform = rayleightest(trial_data.theta_radians) # test expects radians ?


    # return p_uniform
    return None

sub_cond_trial_df = []

for subj_id in av_manifold_df.subj_id.unique():
    for condition in av_manifold_df.condition.unique():

        session_data = av_manifold_df.loc[(av_manifold_df.subj_id == subj_id) & (av_manifold_df.condition == condition)].copy()

        # get all the epoch n trials across epochs and plot those vectors according to r and theta
        # will need to do this via FacetGrid or subplots

        # g=sns.FacetGrid(session_data, col='shifted_epoch_trial')
        # g.map(plot_compass, "r", "theta_radians", "shifted_epoch_trial", r_max=session_data.r.max())



        for shifted_epoch_trial in session_data.shifted_epoch_trial.unique()[1:]: # don't need first one (diff is nan)



            trial_data = session_data.loc[session_data.shifted_epoch_trial == shifted_epoch_trial][['r', 'theta_radians', 'shifted_epoch_trial']].reset_index(drop=True)

            # fig, axes = plt.subplots(2, 2, subplot_kw=dict(polar=True))

            # # should be able to call g.map() using custom plotting fn
            # need to melt to longform first (or refraine from nesting... )
            #
            plot_compass(trial_data.r,trial_data.theta_radians, trial_data.shifted_epoch_trial, r_max=trial_data.r.max())
            #

            # evaluate uniformity of vector direction

            p_uniform, reject_uniform = evaluate_directional_uniformity(trial_data)

            trial_data['p_uniform'] = p_uniform
            trial_data['reject_uniform'] = reject_uniform

            sub_cond_trial_df.append(trial_data)




sub_trial_unif_dfs = []


for subj_id in av_manifold_df.subj_id.unique():


    sub_data = av_manifold_df.loc[(av_manifold_df.subj_id == subj_id)].copy()

    # get all the epoch n trials across epochs and plot those vectors according to r and theta
    # will need to do this via FacetGrid or subplots

    # g=sns.FacetGrid(session_data, col='shifted_epoch_trial')
    # g.map(plot_compass, "r", "theta_radians", "shifted_epoch_trial")



    for shifted_epoch_trial in sub_data.shifted_epoch_trial.unique()[1:]: # don't need first one (diff is nan)



        trial_data = sub_data.loc[sub_data.shifted_epoch_trial == shifted_epoch_trial][['r', 'theta_radians', 'shifted_epoch_trial', 'subj_id']].reset_index(drop=True)

        # fig, axes = plt.subplots(2, 2, subplot_kw=dict(polar=True))

        # # should be able to call g.map() using custom plotting fn
        #
        # plot_compass(trial_data.r,trial_data.theta_radians, r_max=trial_data.r.max())
        #

        # evaluate uniformity of vector direction

        p_uniform, reject_uniform = evaluate_directional_uniformity(trial_data)

        sub_trial_unif_df = pd.DataFrame()

        # add to each plot the value of p

        sub_trial_unif_df['p_uniform'] = p_uniform
        sub_trial_unif_df['reject_uniform'] = reject_uniform
        sub_trial_unif_df['subj_id'] = subj_id
        sub_trial_unif_df['shifted_epoch_trial'] = shifted_epoch_trial


        print(sub_trial_unif_df.head())

        sub_trial_unif_dfs.append(sub_trial_unif_df)

sub_trial_unif_df_all  = pd.concat(sub_trial_unif_dfs, axis=0)



import numpy as np
import pandas as pd
import seaborn as sns
import time
import os
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


av_manifold_df = pd.read_csv('~/Dropbox/loki_0.5/analysis/aggregated_data/av_manifold_df.csv')
av_manifold_df_unit = av_manifold_df.loc[(av_manifold_df.subj_id == 789)].reset_index(drop=True).copy()

n_rows=2; n_cols=5;
n_conditions=9;
plot_mean_components=None
savefig=1
fig_path='/home/krista/Dropbox/loki_0.5/figures/av_manifold_polar_plots'

# from jupyterthemes import jtplot
# jtplot.style(context='talk', fscale=1.1, spines=False, gridlines='--')
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt


def evaluate_directional_uniformity(trial_data, criterion=0.05):
    from astropy.stats import rayleightest

    # plot distribution
    plt.ioff()
    plt.figure()
    plt.title('distribution of angles for epoch trial vector {}:{}'.format((trial_data.shifted_epoch_trial.unique()-1)[0],
     trial_data.shifted_epoch_trial.unique()[0]))
    plt.hist(trial_data.theta_radians)
    plt.xlabel('angle (rad.)'); plt.ylabel('frequency')


    # evaluate statistical sig.
    p_uniform = rayleightest(trial_data.theta_radians) # test expects radians

    reject_uniform = p_uniform <= criterion

    return p_uniform, reject_uniform


# print at max view
# figure out how to check fig size when max.


def plot_compass_v2(av_manifold_df_unit, n_rows=2, ncols=5,
n_conditions=9, plot_mean_components=None, savefig=None, figsize=(10,8), fig_path='/home/krista/Dropbox/loki_0.5/figures/av_manifold_polar_plots'):

    if len(av_manifold_df_unit.condition.unique()) == n_conditions:
        all_conditions = True
    elif len(av_manifold_df_unit.condition.unique()) == 1:
        all_conditions = False
    else:
        raise ValueError('check n_conditions in supplied df')


    fig, main_ax = plt.subplots(n_rows, n_cols, subplot_kw=dict(polar=True), figsize=figsize)
    ax_list = main_ax.flatten()
    kw = dict(arrowstyle="->", color='gray', lw=2)
    mean_kw = dict(arrowstyle="->", color='black', lw=2.5)

    shifted_epoch_trials = np.sort(av_manifold_df_unit.shifted_epoch_trial.unique())[1:]

    for shifted_trial, ax in zip(shifted_epoch_trials, ax_list):
        print(shifted_trial)
        data = av_manifold_df_unit.loc[av_manifold_df_unit.shifted_epoch_trial == shifted_trial].reset_index(drop=True).copy()
        print(data.head())
        r_max = data.r.max()
        print(r_max)
        [ax.annotate("", xy=(theta, rad), xytext=(0,0), arrowprops=kw) for rad, theta in zip(data.r, data.theta_radians)]
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.tick_params(axis='y', which='major', labelsize=8)
        ax.set_ylim(0, r_max)

        mean_rad = np.mean(data.r)
        mean_theta = np.mean(data.theta_radians)

        if plot_mean_components:
            rads = np.arange(0, (2*np.pi), 0.01)
            mean_theta_arr = np.ones_like(rads)*mean_theta
            plt.polar(mean_theta, r_max,  '.', color='black', clip_on=False, markersize=20)
            plt.polar(rads, mean_theta_arr, '-', color='black')
        ax.annotate("", xy=(mean_theta, mean_rad), xytext=(0,0), arrowprops=mean_kw)

        if shifted_trial == 0:
            ax.set_title(r'$\vec{{_{{-\textrm{}:\textrm{}}}}}$'.format(abs(int(shifted_trial-1)), int(shifted_trial)), fontsize=20, y=1.10)
        else:
            ax.set_title(r'$\vec{{_{{\textrm{}:\textrm{}}}}}$'.format(int(shifted_trial)-1, int(shifted_trial)), fontsize=20, y=1.10)

        p_uniform, reject_uniform = evaluate_directional_uniformity(data)

        if reject_uniform:
            ax.annotate('* p = ' + str(np.round(p_uniform, 3)), xy=(0, r_max), xytext=(1.22,r_max+0.01), fontsize=18, color='tomato')

        print('trial {} p_uniform {}'.format(shifted_trial, p_uniform))
        print('reject uniform dist.? {}'.format( reject_uniform))
        # time.sleep(1)


    fig.delaxes(ax_list[-1]) # odd number of plots so delete last subplot
    fig.subplots_adjust(top = 0.90, bottom=0.01, hspace=0.03, wspace=0.2)

    if all_conditions:
        fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + 'all_conditions_av_manifold_polar.png')
        fig.suptitle('subject ' + str(int(av_manifold_df_unit.subj_id.unique()[0])) +  ': all data', fontsize=20)
    else:
        fig_name = (str(int(av_manifold_df_unit.subj_id.unique()[0])) + '_' + str(int(av_manifold_df_unit.condition.unique()[0])) +'_av_manifold_polar.png')
        fig.suptitle('subject ' + str(int(av_manifold_df_unit.subj_id.unique()[0])) + ': ' +
        '$\lambda = $' + str(int(av_manifold_df_unit.lambda_val.unique()[0])) + 'p = ' + int(str(av_manifold_df_unit.p_optimal.unique()[0])), fontsize=20)
        fig.show()

    if savefig:
        # fig.tight_layout()
        # figManager = plt.get_current_fig_manager()
        # figManager.window.showMaximized()
        # fig.show()
        fig.savefig(os.path.join(fig_path, fig_name))

    return fig, main_ax, p_uniform_list, reject_uniform_list
#
# for subj_id in av_manifold_df.subj_id.unique():
#     sub_data = av_manifold_df.loc[av_manifold_df.subj_id == subj_id].reset_index(drop=True).copy()



# _____

# epoch_window = 10
#
# assert len(reg_est_df.epoch_trial.unique()) == epoch_window, 'check epoch len'
#
# for subj_id in reg_est_df.subj_id.unique():
#
#     sub_data = reg_est_df.loc[reg_est_df.subj_id == subj_id]
#
#     x = sub_data.a_est.values
#     y = sub_data.v_est.values
#
#     print(x, y)
#
#     fig, ax = plt.subplots()
#     plt.xticks(rotation=45)
#     plt.title('subject ' + str(subj_id))
#     ax.set_ylabel(r'$\hat{v}$')
#     ax.set_xlabel(r'$\hat{a}$')
#     plt.tight_layout()
#     x_init, y_init = [],[]
#     line, = ax.plot(x_init, y_init, '--', markersize=8)
#
    #
    # annotation = ax.annotate(
    #     '', xy=(1,0), xytext=(-1,0),
    #     arrowprops = {'arrowstyle': "->"},
    # )
    #
    # # # create some random data
    # # n_samples = 5000
    # # base_x = np.random.randn(n_samples)
    # # base_y = base_x
    #
    # # y = np.sin(base_x**2 + base_y**2)
    # # x = np.cos(base_x*base_y)
    #
    # # TODO: maybe just use the annotation object instead of line
    #
    # plt.xlim(x.min()-.0001, x.max()+.0001)
    # plt.ylim(y.min()-.012, y.max()+.012)
    #
    # def init():
    #     line.set_data([], [])
    #     return line,
    #
    # def animate(i):
    #     line.set_data(x[:i],y[:i])
    #     annotation.set_position((x[i-1],y[i-1]))
    #     annotation.xy = (x[i-1],y[i-1])
    #     if i > 0:
    #         annotation.set_text('    t' + str(i-2))
    #
    #     return line, annotation
    #
    # anim = FuncAnimation(fig, animate, init_func=init,
    #                                blit=True,  frames=len(x)+1,
    #                                interval=300, repeat=True) # reduce interval to speed up
    # anim.save(str(subj_id) + '_av_trial.mp4', fps=1, dpi=600, bitrate=-1) # bitrate=-1 lets mpl choose
