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
from spinncer.cerebellum import Cerebellum

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
        results_file += ".npz"

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
    # plot order
    plot_order = ['glomerulus',
                  'granule', 'golgi',
                  'stellate', 'basket',
                  'purkinje', 'dcn']
    n_plots = float(len(plot_order))
    # Plotting results for ...
    print("Plotting results for", results_file)

    # Retrieve information from results file
    all_spikes = data['all_spikes'].ravel()[0]
    final_connectivity = data['final_connectivity'].ravel()[0]
    all_neurons = data['all_neurons'].ravel()[0]
    sim_params = data['simulation_parameters'].ravel()[0]
    simtime = data['simtime'] * ms
    timestep = sim_params['argparser']['timestep'] * ms
    stimulus_params = data['stimulus_params'].ravel()[0]
    starts = np.cumsum(np.concatenate(([0], stimulus_params['stim_times'])))
    time_filter = starts
    stim_durations = sim_params['argparser']['stim_times']
    stimulus_periods = len(stim_durations)
    filtered_firing_rates = {}
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
        # temporary variable to store the population level firing rates
        # before, during and after stimulation
        _filtered_spike_rates = np.zeros(stimulus_periods)
        _spike_times = spikes[:, 1]
        for period in range(stimulus_periods):
            # TODO enhance this to produce histograms of firing rates / pop
            _filtered_spike_times = np.logical_and(
                _spike_times >= time_filter[period],
                _spike_times < time_filter[period + 1])
            _filtered_spike_rates[period] = \
                np.count_nonzero(_filtered_spike_times) / \
                (stim_durations[period] * ms)
        # save the firing rate for the average neuron in this population
        filtered_firing_rates[pop] = _filtered_spike_rates / all_neurons[pop]
    # Report average firing rates before, during and after stimulation
    print("=" * 60)
    print("Average firing rates before, during and after stimulation")
    print("-" * 60)
    for pop, spikes in all_spikes.items():
        _x = filtered_firing_rates[pop] / Hz
        before = _x[0]
        during = _x[1]
        after = _x[2]
        print("\t{:10}->[{:>8.2f}, {:>8.2f}, {:>8.2f}] Hz".format(
            pop, before, during, after), "per neuron")
    print("=" * 60)
    # Report weights values
    print("Average weight per projection")
    print("-" * 60)
    for key in final_connectivity:
        # Connection holder annoyance here:
        conn = np.asarray(final_connectivity[key])
        if conn.size == 0:
            print("Skipping analysing connection", key)
            continue
        if len(conn.shape) == 1 or conn.shape[1] != 4:
            try:
                x = np.concatenate(conn)
                conn = x
            except:
                pass
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
        # replace with percentage of difference
        original_conn = np.abs(CONNECTIVITY_MAP[key]["weight"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95
        _c = Fore.GREEN if is_close else Fore.RED

        print("{:10} -> {}{:4.8f}{} uS".format(
            key, _c, mean, Style.RESET_ALL),
            "c.f. {: 4.8f} uS ({:>7.2%})".format(
                CONNECTIVITY_MAP[key]["weight"], proportion))
    print("=" * 60)
    print("Plotting figures...")
    print("-" * 60)
    # raster plot
    print("Plotting spiking raster plot for each population")
    f, axes = plt.subplots(len(all_spikes.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=700)
    for index, pop in enumerate(plot_order):
        axes[index].scatter(all_spikes[pop][:, 1],
                            all_spikes[pop][:, 0],
                            color=viridis_cmap(index / (n_plots + 1)),
                            s=.5)
        axes[index].set_title(pop)
        # axes[index].set_ylim([0, all_neurons[pop]])
    plt.xlabel("Time (ms)")
    plt.savefig(os.path.join(sim_fig_folder,
                             "raster_plots.png"))
    plt.close(f)

    # plot .1 ms PSTH
    print("Plotting PSTH for each timestep")
    f, axes = plt.subplots(len(spikes_per_timestep.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(plot_order):
        axes[index].bar(np.arange(spikes_per_timestep[pop].size),
                        spikes_per_timestep[pop],
                        color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "timestep_psth.png"))
    plt.close(f)

    # plot sorted .1 ms PSTH
    print("Plotting sorted PSTH for each timestep")
    f, axes = plt.subplots(len(plot_order), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(spikes_per_timestep.keys()):
        axes[index].bar(np.arange(spikes_per_timestep[pop].size),
                        np.sort(spikes_per_timestep[pop]),
                        color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "sorted_timestep_psth.png"))
    plt.close(f)

    # plot 3 ms PSTH
    print("Plotting PSTH in bins of 3 ms")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(plot_order):
        axes[index].bar(np.arange(spikes_per_3ms[pop].size), spikes_per_3ms[pop],
                        color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "timestep_psth_3ms.png"))
    plt.close(f)

    # plot sorted 3 ms PSTH
    print("Plotting sorted PSTH in bins of 3 ms")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(plot_order):
        axes[index].bar(np.arange(spikes_per_3ms[pop].size),
                        np.sort(spikes_per_3ms[pop]),
                        color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "sorted_timestep_psth_3ms.png"))
    plt.close(f)

    # TODO plot weight histogram

    # TODO plot centred connectivity

    # plot as gold standard
    # Select bin size
    n_bins = int(round(simtime / ms / 20))
    print("Plotting to gold standard")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=500)
    for index, pop in enumerate(plot_order):
        axes[index].hist(all_spikes[pop][:, 1], n_bins,
                         color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "gold_standard_psth_3ms.png"))
    plt.close(f)

    print("=" * 60)


if __name__ == "__main__":
    # Constants
    fig_folder = "figures"

    # res = "results/network_results_500x_ssa_55um"
    # spike_analysis(res, fig_folder)

    # res = "results/spikes_by_beatrice"
    # spike_analysis(res, fig_folder)

    # Analyse runs below
    res = "results/gold_standards/gold_standard_results_158"
    spike_analysis(res, fig_folder)

    res = "results/gold_standards/gold_standard_results_400"
    spike_analysis(res, fig_folder)
