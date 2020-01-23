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
from spinncer.utilities.constants import CONNECTIVITY_MAP
from spinncer.utilities.neo_converter import convert_spikes
from colorama import Fore, Style, init as color_init

mlib.use('Agg')
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
    other_recordings = data['other_recordings'].ravel()[0]
    simtime = data['simtime'] * ms
    timestep = sim_params['argparser']['timestep'] * ms
    stimulus_params = data['stimulus_params'].ravel()[0]
    starts = np.cumsum(np.concatenate(([0], stimulus_params['stim_times'])))
    time_filter = starts
    stim_durations = sim_params['argparser']['stim_times']
    stimulus_periods = len(stim_durations)
    average_firing_rates = {}
    filtered_firing_rates = {}
    # Check if using neo blocks
    neo_all_spikes = {}
    for pop, potential_neo_block in all_spikes.items():
        if isinstance(potential_neo_block, neo.Block):
            # make a copy of the spikes dict
            neo_all_spikes[pop] = potential_neo_block
            all_spikes[pop] = convert_spikes(potential_neo_block)
    # Report useful parameters
    print("=" * 80)
    print("Simulation parameters")
    print("-" * 80)
    pp(sim_params)
    # Report useful parameters
    print("=" * 80)
    print("Analysis report")
    print("-" * 80)
    print("Current time",
          plt.datetime.datetime.now().strftime("%H:%M:%S on %d.%m.%Y"))
    # Report number of neurons
    print("=" * 80)
    print("Number of neurons in each population")
    print("-" * 80)
    for pop in plot_order:
        print("\t{:10} -> {:10} neurons".format(pop, all_neurons[pop]))
    # Pre-compute conversions
    time_to_bin_conversion = 1. / (timestep / ms)
    bins_in_3ms = int(3 * time_to_bin_conversion)
    no_timesteps = int(simtime / ms * time_to_bin_conversion)
    pad_to_compute_3ms_bins = bins_in_3ms - no_timesteps % bins_in_3ms
    # Bincount
    spikes_per_timestep = {}
    spikes_per_3ms = {}
    # Per neuron firing rate in each stimulus period
    per_neuron_firing = {}
    per_neuron_spike_count = {}
    print("=" * 80)
    print("Maximum number of generated spikes per timestep")
    print("-" * 80)
    for pop in plot_order:
        spikes = all_spikes[pop]

        spikes_per_timestep[pop] = \
            np.bincount((spikes[:, 1] * time_to_bin_conversion).astype(int),
                        minlength=no_timesteps)
        max_spikes = np.max(spikes_per_timestep[pop])
        print("\t{:10}->{:6} = {:1.4f}".format(pop, max_spikes,
                                               max_spikes / all_neurons[pop]),
              "per neuron")
        padded_bincount = np.pad(
            spikes_per_timestep[pop],
            (0, pad_to_compute_3ms_bins -
             (spikes_per_timestep[pop].size - no_timesteps)),  # This should be 0 or 1 corresponding to SpiNNaker or NEST
            'constant', constant_values=0)

        reshaped_bincount = padded_bincount.reshape(
            int(padded_bincount.shape[0] / bins_in_3ms), bins_in_3ms)

        spikes_per_3ms[pop] = np.sum(reshaped_bincount, axis=1)
        # temporary variable to store the population level firing rates
        # before, during and after stimulation
        _filtered_spike_rates = np.zeros(stimulus_periods)
        _spike_times = spikes[:, 1]
        # Initialise per_neuron_firing
        per_neuron_firing[pop] = np.ones((all_neurons[pop],
                                          stimulus_periods)) * -10
        per_neuron_spike_count[pop] = np.ones((all_neurons[pop],
                                               stimulus_periods)) * -10
        for period in range(stimulus_periods):
            _filtered_spike_times = np.logical_and(
                _spike_times >= time_filter[period],
                _spike_times < time_filter[period + 1])
            _filtered_spike_rates[period] = \
                np.count_nonzero(_filtered_spike_times) / \
                (stim_durations[period] * ms)
            for nid in range(all_neurons[pop]):
                _spikes_for_nid = spikes[spikes[:, 0] == nid][:, 1]
                _no_spike_for_nid = np.count_nonzero(np.logical_and(
                    _spikes_for_nid >= time_filter[period],
                    _spikes_for_nid < time_filter[period + 1]))
                per_neuron_spike_count[pop][nid, period] = _no_spike_for_nid
                per_neuron_firing[pop][nid, period] = \
                    _no_spike_for_nid / (stim_durations[period] * ms)
        # save the firing rate for the average neuron in this population
        average_firing_rates[pop] = _filtered_spike_rates / all_neurons[pop]
    # Report average firing rates before, during and after stimulation
    print("=" * 80)
    print("Average firing rates before, during and after stimulation")
    print("-" * 80)
    for pop in plot_order:
        _x = average_firing_rates[pop] / Hz
        before = _x[0]
        during = _x[1]
        after = _x[2]
        print("\t{:10}->[{:>8.2f}, {:>8.2f}, {:>8.2f}] Hz".format(
            pop, before, during, after), "per neuron")
    print("=" * 80)

    print("=" * 80)
    print("FILTERED average firing rates before, during and after stimulation")
    print("-" * 80)
    for pop in plot_order:
        # Retrieve firing rate for each neuron in each of the periods
        _x = per_neuron_firing[pop]
        # _filtered_rates = np.zeros(stimulus_periods)
        _excited_map = np.zeros(_x.shape[0])
        _inhibited_map = np.zeros(_x.shape[0])

        # filter out neurons that are not "excited" as defined in the
        # scaffold paper, i.e., cells that fire twice as much when stimulated as
        # compared to before stimulation
        _excited_map = np.greater(_x[:, 1], 2 * _x[:, 0])
        # See SC email from 23/01/2020
        # baseline firing rate b
        # stimulation firing rate s
        # correct: b > 2 * s
        _inhibited_map = np.greater(_x[:, 0], 2 * _x[:, 1])

        # filter out neurons that don't fire at all
        _neurons_that_fire = per_neuron_spike_count[pop][:, 1] > 0
        _excited_map = np.logical_and(
            _excited_map, _neurons_that_fire
        )
        _inhibited_map = np.logical_and(
            _inhibited_map, _neurons_that_fire
        )

        if pop == "granule":
            # for GrC filter out neurons that don't fire more than once
            _excited_map = np.logical_and(
                _excited_map, per_neuron_spike_count[pop][:, 1] > 1
            )

        excited_filteread_mean = np.mean(_x[_excited_map], axis=0)
        excited_filteread_std = np.std(_x[_excited_map], axis=0)
        excited_before = excited_filteread_mean[0]
        excited_during = excited_filteread_mean[1]
        excited_after = excited_filteread_mean[2]
        excited_before_std = excited_filteread_std[0]
        excited_during_std = excited_filteread_std[1]
        excited_after_std = excited_filteread_std[2]
        print("\t{:10} excited   ->[{:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}] Hz".format(
            pop,
            excited_before, excited_before_std,
            excited_during, excited_during_std,
            excited_after, excited_after_std),
            "per neuron")
        no_excited = np.count_nonzero(_excited_map)
        print("\t\t\t {:6d} excited neurons, i.e. {:7.2%} of cells".format(
            no_excited, no_excited / all_neurons[pop]))

        inhibited_filteread_mean = np.mean(_x[_inhibited_map], axis=0)
        inhibited_filteread_std = np.std(_x[_inhibited_map], axis=0)
        inhibited_before = inhibited_filteread_mean[0]
        inhibited_during = inhibited_filteread_mean[1]
        inhibited_after = inhibited_filteread_mean[2]
        inhibited_before_std = inhibited_filteread_std[0]
        inhibited_during_std = inhibited_filteread_std[1]
        inhibited_after_std = inhibited_filteread_std[2]
        print("\t{:10} inhibited ->[{:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}] Hz".format(
            pop,
            inhibited_before, inhibited_before_std,
            inhibited_during, inhibited_during_std,
            inhibited_after, inhibited_after_std),
            "per neuron")
        no_inhibited = np.count_nonzero(_inhibited_map)
        print("\t\t\t {:6d} inhibited neurons, i.e. {:7.2%} of cells".format(
            no_inhibited, no_inhibited / all_neurons[pop]))
        print("\t{:10} neurons didn't fire at all".format(
            np.count_nonzero(np.invert(_neurons_that_fire))))
        print("-" * 80)
    print("=" * 80)

    if analysis_args.worst_case_spikes:
        # Count incoming spikes only if we care -- this takes a while
        inc_spike_count = {k: np.zeros((all_neurons[k], no_timesteps + 1)) for k in all_neurons.keys()}

    # flag set if some connectivity exists
    conn_exists = False

    # Report weights values
    print("Average weight per projection")
    print("-" * 80)
    conn_dict = {}
    for key in final_connectivity:
        # Connection holder annoyance here:
        conn = np.asarray(final_connectivity[key])
        if final_connectivity[key] is None or conn.size == 0:
            print("Skipping analysing connection", key)
            continue
        conn_exists = True
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
        conn_dict[key] = conn
        mean = np.mean(conn[:, 2])
        # replace with percentage of difference
        original_conn = np.abs(CONNECTIVITY_MAP[key]["weight"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        # assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95
        _c = Fore.GREEN if is_close else Fore.RED

        print("{:10} -> {}{:4.8f}{} uS".format(
            key, _c, mean, Style.RESET_ALL),
            "c.f. {: 4.8f} uS ({:>7.2%})".format(
                CONNECTIVITY_MAP[key]["weight"], proportion))

    # Check voltage information

    print("=" * 80)
    print("Input current analysis")
    print("-" * 80)
    all_voltages = {}

    if other_recordings is not None:
        for pop in plot_order:
            try:
                curr_v = other_recordings[pop]
                if curr_v is None:
                    raise KeyError()
                else:
                    curr_v = curr_v['v']
            except KeyError:
                print("No voltage information for", pop)
                continue
            # Create a useful, aligned numpy array of shape (sh0, sh1)
            # where sh0 = number of neurons in pop
            # and   sh1 = number of timesteps in the simulation
            if isinstance(curr_v, neo.Block):
                try:
                    all_voltages[pop] = np.array(curr_v.segments[0].filter(name='v')[0]).T
                    if pop == "dcn":
                        # TODO plot voltage traces
                        pass
                except AttributeError as ae:
                    print("[WARNING]", pop, "has no voltage information.")
                    all_voltages[pop] = np.ones((1, no_timesteps)) * 3.14

            elif curr_v.shape[1] == 4:
                all_voltages[pop + "_exc"] = np.zeros((all_neurons[pop], no_timesteps))
                all_voltages[pop + "_inh"] = np.zeros((all_neurons[pop], no_timesteps))
                for nid, time, v_exc, v_inh in curr_v:
                    all_voltages[pop + "_exc"][int(nid), int(time * time_to_bin_conversion)] = v_exc
                    all_voltages[pop + "_inh"][int(nid), int(time * time_to_bin_conversion)] = v_inh
            else:
                print("Only synaptic current contribution is supported "
                      "for analysis currently")
        # Report statistics here
        for key, v in all_voltages.items():
            nid, timestep = np.unravel_index(np.argmax(v, axis=None), v.shape)
            print("{:10}-> neuron {:>8d} received {:>6d}".format(
                key, int(nid), int(np.max(v))),
                "nA in timestep #{:8d}".format(int(timestep)))
    else:
        print("No voltage information present.")

    if analysis_args.worst_case_spikes:
        # The following is expensive time wise
        for key, conn in conn_dict.items():
            post_pop = CONNECTIVITY_MAP[key]["post"]
            pre_pop = CONNECTIVITY_MAP[key]["pre"]
            curr_spikes = all_spikes[pre_pop]
            for nid, t in curr_spikes:
                nid = int(nid)
                times = int(t * time_to_bin_conversion)
                targets = conn[conn[:, 0] == nid][:, 1].astype(int)
                inc_spike_count[post_pop][targets, times] += 1

    if conn_exists and analysis_args.worst_case_spikes:
        print("=" * 80)
        print("Incoming spikes statistics")
        print("-" * 80)
        for pop, counts in inc_spike_count.items():
            nid, timestep = np.unravel_index(np.argmax(counts, axis=None), counts.shape)
            print("{:10}-> neuron {:>8d} saw {:>6d}".format(
                pop, int(nid), int(np.max(counts))),
                "spikes in timestep #{:8d}".format(int(timestep)))
            # Break down spikes in that timestep

    print("=" * 80)
    print("Plotting figures...")
    print("-" * 80)
    # plot firing rate histogram per PSTH region
    print("Plotting firing rate histograms")
    f, axes = plt.subplots(len(plot_order), 3,
                           figsize=(14, 20), sharex=True, dpi=700)
    for index, pop in enumerate(plot_order):
        for period in range(stimulus_periods):
            curr_ax = axes[index, period]
            curr_ax.hist(per_neuron_firing[pop][:, period],
                         color=viridis_cmap(index / (n_plots + 1)),
                         bins=20)
            if period == 1:
                curr_ax.set_title(pop)
            curr_ax.set_xlabel("Hz")
            curr_ax.xaxis.set_tick_params(which='both',
                                          labelbottom=True)
            curr_ax.set_xticks([0, 75, 150])

    plt.savefig(os.path.join(sim_fig_folder,
                             "neuron_firing_rate_hist.png"))
    plt.close(f)

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
                           figsize=(14, 20), sharex=True, dpi=700)
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
                           figsize=(14, 20), sharex=True, dpi=700)
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
                           figsize=(14, 20), sharex=True, dpi=700)
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
                           figsize=(14, 20), sharex=True, dpi=700)
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
                           figsize=(14, 20), sharex=True, dpi=700)
    for index, pop in enumerate(plot_order):
        axes[index].hist(all_spikes[pop][:, 1], n_bins,
                         color=viridis_cmap(index / (n_plots + 1)))
        axes[index].set_title(pop)
    plt.savefig(os.path.join(sim_fig_folder,
                             "gold_standard_psth_3ms.png"))
    plt.close(f)

    print("=" * 80)


if __name__ == "__main__":
    from spinncer.analysis_argparser import *

    if analysis_args.input and len(analysis_args.input) > 0:
        for in_file in analysis_args.input:
            spike_analysis(in_file, analysis_args.figures_dir)
    else:
        # Constants
        fig_folder = "figures"

        res = "results/gold_standards/gold_standard_results_400_stim_radius_140"
        spike_analysis(res, fig_folder)

        res = "results/gold_standards/gold_standard_results_400_stim_radius_70"
        spike_analysis(res, fig_folder)
