import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm as cm_mlib
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy
from matplotlib import animation, rc, colors
from brian2.units import *
import matplotlib as mlib
from scipy import stats
from pprint import pprint as pp
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
import traceback
import os
import copy
import neo
from datetime import datetime
import warnings
import ntpath
from spinncer.utilities.constants import CONNECTIVITY_MAP, CELL_PARAMS, SPIKE_SOURCE_NAMES
from spinncer.utilities.neo_converter import convert_spikes
from colorama import Fore, Style, init as color_init
import pandas as pd
import string
from matplotlib.ticker import MultipleLocator
import elephant
import quantities as pq

mlib.use('Agg')
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ensure we use viridis as the default cmap
plt.viridis()

# ensure we use the same rc parameters for all matplotlib outputs
mlib.rcParams.update({'font.size': 24})
mlib.rcParams.update({'errorbar.capsize': 5})
mlib.rcParams.update({'figure.autolayout': True})
viridis_cmap = mlib.cm.get_cmap('viridis')

PREFFERED_ORDER = ['mossy_fibers',
                   'glomerulus',
                   'granule', 'golgi',
                   'stellate', 'basket',
                   'purkinje', 'dcn']


def color_for_index(index, size, cmap=viridis_cmap):
    return cmap(index / (size + 1))


def write_sep():
    print("=" * 80)


def write_line():
    print("-" * 80)


def write_header(msg):
    write_sep()
    print(msg)
    write_line()


def write_short_msg(msg, value):
    print("{:40}:{:39}".format(msg, str(value)))


def write_value(msg, value):
    print("{:60}:{:19}".format(msg, str(value)))


def get_plot_order(for_keys):
    # Compute plot order
    plot_order = []
    # only focus on keys for pops that have spikes
    key_duplicate = list(for_keys)
    key_duplicate.sort()
    for pref in PREFFERED_ORDER:
        for i, key in enumerate(key_duplicate):
            if pref in key:
                plot_order.append(key)
                key_duplicate.pop(i)

    # add remaining keys to the end of plot order
    plot_order += key_duplicate
    return plot_order


def scatter_hist(x, y, ax, ax_histx, ax_histy):
    """
    https://matplotlib.org/gallery/lines_bars_and_markers/scatter_hist.html
    """
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')


def imshow_hist(x, ax, ax_histx, ax_histy, only_count_1s=False):
    """
        based on https://matplotlib.org/gallery/lines_bars_and_markers/scatter_hist.html
    @param x:
    @param ax:
    @param ax_histx:
    @param ax_histy:
    @param only_count_1s:
    @return:
    """
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ax.imshow(x)
    if only_count_1s:
        x[x != 1] = 0

    # now determine nice limits by hand:
    binwidth = 1
    xymax = np.nanmax(np.nanmax(np.abs(x)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    x_dim = x.shape[0]
    y_dim = x.shape[1]

    ax_histx.bar(np.arange(y_dim),
                 height=np.nansum(x, axis=0),
                 width=1.,
                 )

    ax_histy.barh(np.arange(x_dim),
                  height=1.,
                  width=np.nansum(x, axis=1),
                  )
    ax.set_xlim([0, x.shape[1]])
    ax.set_ylim([0, x.shape[0]])


COMMON_DISPLAY_NAMES = {
    'f_peak': "$f_{peak}$ (Hz)",
    'spinnaker': "SpiNNaker",
    'nest': "NEST",
    'stim_radius': "Stimulation radius ($\mu m$)",
    'glomerulus cells': "Glom",
    'glomerulus': "Glom",
    'mossy_fibers': "MF",
    'mossy_fibres': "MF",
    'climbing_fibres': "CF",
    'climbing_fibers': "CF",
    'granule cells': "GrC",
    'granule': "GrC",
    'dcn cells': "DCNC",
    'dcn': "DCNC",
    'golgi cells': "GoC",
    'golgi': "GoC",
    'purkinje cells': "PC",
    'purkinje': "PC",
    'stellate cells': "SC",
    'stellate': "SC",
    'basket cells': "BC",
    'basket': "BC",
    'max_spikes_in_a_tick': "Peak # of MC packets",
    'dumped_from_a_link': "Dropped from link",
    'send_multicast_packets': "MC packets sent",
    'max_dmas_in_a_tick': "Peak # of DMAs",
    'max_pipeline_restarts': "Peak # of pipeline restarts",
    'router_provenance': "Router",
    # Connection names
    'aa_goc': "aa-GoC",
    'aa_pc': "aa-PC",
    'bc_pc': "BC-PC",
    'gj_bc': "BC-BC",
    'gj_goc': "GoC-GoC",
    'gj_sc': "SC-SC",
    'glom_dcn': "Glom-DCNC",
    'glom_goc': "Glom-GoC",
    'glom_grc': "Glom-GrC",
    'goc_grc': "GoC-GrC",
    'pc_dcn': "PC-DCNC",
    'pf_bc': "pf-BC",
    'pf_goc': "pf-GoC",
    'pf_pc': "pf-PC",
    'pf_sc': "pf-SC",
    'sc_pc': "SC-PC",
    'mf_glom': "MF-Glom"
}


def capitalise(name):
    return string.capwords(
        " ".join(
            str.split(name, "_")
        ))


def use_display_name(name):
    name = name.lower()
    return COMMON_DISPLAY_NAMES[name] \
        if name in COMMON_DISPLAY_NAMES.keys() \
        else capitalise(name)


def save_figure(plt, name, extensions=(".png",), **kwargs):
    for ext in extensions:
        write_short_msg("Plotting", name + ext)
        plt.savefig(name + ext, **kwargs)


def is_iterable(obj):
    """
    https://stackoverflow.com/questions/1952464/in-python-how-do-i-determine-if-an-object-is-iterable
    """
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True


# Composable functions to easily get the necessary spiking information
def filter_spikes_by_post_ids(pre_spikes, connection, post_ids):
    if not is_iterable(post_ids):
        post_ids = np.array([post_ids]).astype(int)
    pre_ids = connection[np.in1d(connection[:, 1].astype(int), post_ids)][:, 0].astype(int)
    return pre_spikes[np.in1d(pre_spikes[:, 0].astype(int), pre_ids)]


def filter_spikes_by_timestep(spikes, timestep, dt=0.1):
    """
    Call filter_spikes_by_times
    """
    return filter_spikes_by_time(spikes, timestep * dt)


def filter_spikes_by_time(spikes, time):
    return spikes[np.isclose(spikes[:, 1], time)]



# TODO make this self-contained
def plot_3d_cell_positions(position, plot_order,
                           highlight_afferents_for=None, conn_params=None, final_connectivity=None,
                           hightlight_rasterisation_for=None, az=0, rot=0, scenic=False):
    # Do ID mapping for pre and post in highlight_afferents_for here
    h_pre_ids = []
    h_post_ids = []
    h_active_pre_ids = []
    if highlight_afferents_for:
        h_post_ids = highlight_afferents_for['post_ids']
        timestep = highlight_afferents_for['timestep']
        curr_proj = highlight_afferents_for['proj']
        curr_conn = final_connectivity[highlight_afferents_for['proj']]
        for hpi in h_post_ids:
            # Also retrieve active pre ids
            # This involves computing active afferents
            passive_pres = curr_conn[curr_conn[:, 1]==hpi][:, 0].astype(int)
            if timestep:
                active_pres = filter_spikes_by_timestep(all_spikes[conn_params[curr_proj]['pre']], timestep)
                # highligh only active cells
                h_pre_ids.extend(np.intersect1d(passive_pres, active_pres))
            else:
                h_pre_ids.extend(passive_pres)


    h_pre_ids = np.asarray(h_pre_ids).ravel()

    # Order position ids based on plot_order
    id_order = [POPULATION_ID[po] for po in plot_order]
    #
    fig = plt.figure(figsize=(10,8), dpi=300)
    ax = fig.add_subplot(111, projection='3d')

    for ind, (_id, pop) in enumerate(zip(id_order, plot_order)):
        if selective_display and not to_display[pop]:
            continue
        curr_data = position[position['pop']==_id]
        curr_alpha= alphas[pop]
        curr_marker = markers[pop]
        if highlight_afferents_for and pop == conn_params[highlight_afferents_for['proj']]['post']:
            print("Highlight afferents for", pop, highlight_afferents_for)
            # Reduced alpha for non-selected ids
            selected_points = curr_data.iloc[h_post_ids]
            ax.scatter(selected_points.x, selected_points.y, selected_points.z,
                       color=viridis_cmap(ind / (n_plots + 1)),
                       alpha=1, marker="x", label="Selected " + use_display_name(pop))

            if pop != conn_params[highlight_afferents_for['proj']]['pre']:
                unselected_points = curr_data.loc[curr_data.index.difference(h_post_ids)]
                ax.scatter(unselected_points.x, unselected_points.y, unselected_points.z,
                           color=viridis_cmap(ind / (n_plots + 1)),
                           alpha=low_alpha, marker=curr_marker,
                           label="Unselected " + use_display_name(pop))

        if highlight_afferents_for and pop == conn_params[highlight_afferents_for['proj']]['pre']:
            print("Highlight afferents for", pop, highlight_afferents_for)
            # need to run scatter plot twice. Once for selected. Once for ignored.

            selected_points = curr_data.iloc[h_pre_ids]
            ax.scatter(selected_points.x, selected_points.y, selected_points.z,
                       color=viridis_cmap(ind / (n_plots + 1)),
                       alpha=.2, marker=curr_marker,
                       label="Selected " + use_display_name(pop))

            if pop != conn_params[highlight_afferents_for['proj']]['post']:
                unselected_points = curr_data.loc[curr_data.index.difference(h_pre_ids)]
                ax.scatter(unselected_points.x, unselected_points.y, unselected_points.z,
                           color=viridis_cmap(ind / (n_plots + 1)),
                           alpha=low_alpha, marker=curr_marker,
                           label="Unselected " + use_display_name(pop))

        if hightlight_rasterisation_for and pop in hightlight_rasterisation_for:
            id_range = np.asarray(curr_data['nid'])- nid_offset[pop]
            gen_colours = list(map(viridis_cmap, id_range/(all_neurons[pop] + 1)))
            ax.scatter(curr_data.x, curr_data.y, curr_data.z,
                       color=gen_colours,
                       alpha=.8, marker=curr_marker,
                       label=use_display_name(pop))
        elif not highlight_afferents_for:
            ax.scatter(curr_data.x, curr_data.y, curr_data.z,
                       color=viridis_cmap(ind / (n_plots + 1)),
                       alpha=curr_alpha, marker=curr_marker,
                       label=use_display_name(pop))
        # After plotting all of the static connectivity, highlight active connectivity
        if highlight_afferents_for and highlight_afferents_for['timestep']:
            # TODO retrieve ids of active pre connections
            pass

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-5, 405)
    ax.set_ylim(-5, 405)
    # plt.legend(loc="best")
    if scenic:
        ax.set_axis_off()
        plt.tight_layout()
    else:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.view_init(az, rot)  # first number controls azymuth, second control rotation
    save_figure(plt, os.path.join(per_neuron_figs,
                                         "3D_positions_{}_{}".format(az, rot)),
                        extensions=['.png', ])
    plt.show()
    plt.close(fig)
