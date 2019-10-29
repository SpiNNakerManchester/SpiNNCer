import numpy as np
import matplotlib.pyplot as plt
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
from spinncer.utilities.constants import CONNECTIVITY_MAP
from colorama import Fore, Style, init as color_init

warnings.filterwarnings("ignore", category=UserWarning)

# ensure we use viridis as the default cmap
plt.viridis()

# ensure we use the same rc parameters for all matplotlib outputs
mlib.rcParams.update({'font.size': 24})
mlib.rcParams.update({'errorbar.capsize': 5})
mlib.rcParams.update({'figure.autolayout': True})
viridis_cmap = mlib.cm.get_cmap('viridis')


def spike_analysis(results_file, fig_folder):
    # Retrieve results file
    try:
        data = np.load(results_file, allow_pickle=True)
    except FileNotFoundError:
        data = np.load(results_file + ".npz", allow_pickle=True)

    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)

    # Create figures folder for this results_file
    sim_fig_folder = os.path.join(fig_folder,
                                  str(ntpath.basename(results_file))[:-4])
    if not os.path.isdir(sim_fig_folder) and not os.path.exists(sim_fig_folder):
        os.mkdir(sim_fig_folder)
    # Set up colours
    color_init(strip=False)

    # Retrieve information from results file
    all_spikes = data['all_spikes'].ravel()[0]
    final_connectivity = data['final_connectivity'].ravel()[0]
    all_neurons = data['all_neurons'].ravel()[0]
    sim_params = data['simulation_parameters'].ravel()[0]
    simtime = data['simtime'] * ms
    timestep = sim_params['argparser']['timestep'] * ms
    # Report useful parameters
    print("=" * 60)
    print("Simulation parameters")
    print("-" * 60)
    pp(sim_params)
    print("=" * 60)
    # Pre-compute conversions
    time_to_bin_conversion = 1. / (timestep / ms)
    bins_in_3ms = int(3 * time_to_bin_conversion)
    no_timesteps = int(simtime / ms * time_to_bin_conversion)
    pad_to_compute_3ms_bins = bins_in_3ms - no_timesteps % bins_in_3ms
    # Bincount
    spikes_per_timestep = {}
    spikes_per_3ms = {}
    print("Maximum number of generated spikes per timestep")
    print("-" * 60)
    for pop, spikes in all_spikes.items():
        spikes_per_timestep[pop] = \
            np.bincount((spikes[:, 1] * time_to_bin_conversion).astype(int),
                        minlength=no_timesteps)
        max_spikes = np.max(spikes_per_timestep[pop])
        print("\t{:10}->{:6} = {:1.4f}".format(pop, max_spikes,
                                               max_spikes / all_neurons[pop]),
              "per neuron")
        padded_bincount = np.pad(spikes_per_timestep[pop],
                                 (0, pad_to_compute_3ms_bins),
                                 'constant', constant_values=0)
        reshaped_bincount = padded_bincount.reshape(
            int(padded_bincount.shape[0] / bins_in_3ms), bins_in_3ms)

        spikes_per_3ms[pop] = np.sum(reshaped_bincount, axis=1)
    # Report weights values
    print("=" * 60)
    print("Average weight per projection")
    print("-" * 60)
    for key in final_connectivity:
        # Connection holder annoyance here:
        conn = final_connectivity[key]
        if len(conn.shape) == 1 or conn.shape[1] != 4:
            conn = np.concatenate(conn)
            names = [('source', 'int_'),
                     ('target', 'int_'),
                     ('weight', 'float_'),
                     ('delay', 'float_')]
            useful_conn = np.zeros((conn.shape[0], 4), dtype=np.float)
            for i, (n, _) in enumerate(names):
                useful_conn[:, i] = conn[n].astype(np.float)
            final_connectivity[key] = useful_conn.astype(np.float)
            conn = useful_conn.astype(np.float)

        mean = np.mean(conn[:, 2])
        is_close = np.isclose(mean, CONNECTIVITY_MAP[key]["weight"], 1.0e-3)
        _c = Fore.GREEN if is_close else Fore.RED

        print("\t{:10} -> {}{:4.8f}{} uS".format(
            key, _c, mean, Style.RESET_ALL),
            "c.f. {:4.8f} uS".format(CONNECTIVITY_MAP[key]["weight"]))
    print("=" * 60)
    print("Plotting figures...")
    # plot .1 ms PSTH
    print("Plotting PSTH for each timestep")
    f, axes = plt.subplots(len(spikes_per_timestep.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(spikes_per_timestep.keys()):
        axes[index].bar(np.arange(spikes_per_timestep[pop].size),
                        spikes_per_timestep[pop])
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "timestep_psth.png"))
    plt.close(f)

    # plot sorted .1 ms PSTH
    print("Plotting sorted PSTH for each timestep")
    f, axes = plt.subplots(len(spikes_per_timestep.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(spikes_per_timestep.keys()):
        axes[index].bar(np.arange(spikes_per_timestep[pop].size),
                        np.sort(spikes_per_timestep[pop]))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "sorted_timestep_psth.png"))
    plt.close(f)

    # plot 3 ms PSTH
    print("Plotting PSTH in bins of 3 ms")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(spikes_per_3ms.keys()):
        axes[index].bar(np.arange(spikes_per_3ms[pop].size), spikes_per_3ms[pop])
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "timestep_psth_3ms.png"))
    plt.close(f)

    # plot sorted 3 ms PSTH
    print("Plotting sorted PSTH in bins of 3 ms")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(spikes_per_3ms.keys()):
        axes[index].bar(np.arange(spikes_per_3ms[pop].size),
                        np.sort(spikes_per_3ms[pop]))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "sorted_timestep_psth_3ms.png"))
    plt.close(f)

    # TODO plot weight histogram
    print("=" * 60)

    # TODO plot centred connectivity


if __name__ == "__main__":
    # Constants
    fig_folder = "figures"

    # Analyse runs below
    res = "results/spinncer_experiment_130640_29102019"
    spike_analysis(res, fig_folder)

    res = "results/spinncer_no_proj"
    spike_analysis(res, fig_folder)
