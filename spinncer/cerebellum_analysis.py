from spinncer.analysis_common import *

# The following LUT contains the delays for both excitation and inhibition to
# reach that particular population in the case of DCN, or only excitation for
# everything else
DELAY_IN_EXCITATION = {
    'glomerulus': 0,
    'mossy_fibers': 0,
    'granule': 4,
    'granule_cell': 4,
    'dcn': max(4, 10),
    'dcn_cell': max(4, 10),
    'dcn_interneuron': 4,
    'golgi': 4,
    'golgi_cell': 4,
    'purkinje': 6,
    'purkinje_cell': 6,
    'stellate': 9,
    'stellate_cell': 9,
    'basket': 9,
    'basket_cell': 9,
    'io_cell': 0,
}


def color_for_index(index, size, cmap=viridis_cmap):
    return cmap(1 / (size - index + 1))


def plot_analog_signal(data, variable_name, ylabel, plot_order,
                       wanted_times, time_to_bin_conversion,
                       fig_folder,
                       highlight_stim, common_highlight_values):
    print("Plotting {} traces for each population".format(variable_name))
    for index, pop in enumerate(plot_order):
        if pop == "glomerulus":
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            plt.close(f)
            continue
        try:
            values_for_pop = data[pop]
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            # if highlight_stim:
            #     ax = plt.gca()
            #     highlight_area(ax, **common_highlight_values)
            for _ind, _trace in enumerate(values_for_pop):
                plt.plot(_trace,
                         color=color_for_index(_ind, values_for_pop.shape[0]),
                         rasterized=True)
            plt.xticks(wanted_times * time_to_bin_conversion, wanted_times)
            plt.xlabel("Time (ms)")
            plt.ylabel(ylabel)

            save_figure(plt, os.path.join(fig_folder,
                                          "{}_{}".format(pop, variable_name)),
                        extensions=['.png', ])
            plt.close(f)
        except:
            pass


def highlight_area(ax, pop, start, stop, increment):
    _highlight_times = np.arange(start[pop],
                                 stop[pop],
                                 increment)
    ax.fill_between(
        _highlight_times, 0, 1,
        color='grey', alpha=0.3,
        transform=ax.get_xaxis_transform())


def spike_analysis(results_file, fig_folder,
                   worst_case=True, delay_sensitive=False,
                   dark_background=False,
                   highlight_stim=False):
    if dark_background:
        plt.style.use('dark_background')
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

    # Plotting results for ...
    print("=" * 80)
    print("Plotting results for", results_file)
    print("-" * 80)

    # Retrieve information from results file
    all_spikes = data['all_spikes'].ravel()[0]
    try:
        final_connectivity = data['final_connectivity'].ravel()[0]
    except:
        final_connectivity = []
        traceback.print_exc()
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

    conn_params = data['conn_params'].ravel()[0] if 'conn_params' in data.files else CONNECTIVITY_MAP
    cell_params = data['cell_params'].ravel()[0] if 'cell_params' in data.files else CELL_PARAMS
    # Compute plot order
    plot_order = get_plot_order(all_spikes.keys())
    n_plots = float(len(plot_order))
    simulator = sim_params['argparser']['simulator']
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
    pp(cell_params)
    pp(conn_params)
    # Report useful parameters
    print("=" * 80)
    print("Analysis report")
    print("-" * 80)
    print("Current time",
          plt.datetime.datetime.now().strftime("%H:%M:%S on %d.%m.%Y"))
    print("This analysis includes:")
    print("\tUNFILTERED mean firing rates per period (before, during, after "
          "stimulation)")
    print("\tFILTERED   mean firing rates per period (before, during, after "
          "stimulation)")
    print("\traster plots")
    print("\t3ms-binned PSTH")
    print("\ttimestep-binned PSTH")
    if delay_sensitive:
        print("\tfiring rates per period taking into account DELAYS")
    if worst_case:
        print("\tcounting WORST CASE number of afferent spikes per cell")
    # Report number of neurons
    print("=" * 80)
    print("Number of neurons in each population")
    print("-" * 80)
    for pop in plot_order:
        print("\t{:20} -> {:10} neurons".format(pop, all_neurons[pop]))
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

    stim_period_start = {}
    stim_period_end = {}
    per_pop_stim_durations = {k: [] for k in plot_order}
    for pop in plot_order:
        spikes = all_spikes[pop]

        spikes_per_timestep[pop] = \
            np.bincount((spikes[:, 1] * time_to_bin_conversion).astype(int),
                        minlength=no_timesteps)
        max_spikes = np.max(spikes_per_timestep[pop])
        print("\t{:20}->{:6} = {:1.4f}".format(pop, max_spikes,
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

            if delay_sensitive and period == 0:
                time_filter_pre = time_filter[period]
                time_filter_post = time_filter[period + 1] + DELAY_IN_EXCITATION[pop]
            elif delay_sensitive and period == 1:
                time_filter_pre = time_filter[period] + DELAY_IN_EXCITATION[pop]
                time_filter_post = time_filter[period + 1] + DELAY_IN_EXCITATION[pop]
            elif delay_sensitive and period == 2:
                time_filter_pre = time_filter[period] + DELAY_IN_EXCITATION[pop]
                time_filter_post = time_filter[period + 1]
            else:
                time_filter_pre = time_filter[period]
                time_filter_post = time_filter[period + 1]

            if period == 1:
                stim_period_start[pop] = time_filter_pre
                stim_period_end[pop] = time_filter_post
            per_pop_stim_durations[pop].append(time_filter_post - time_filter_pre)
            current_period_duration = per_pop_stim_durations[pop][period]
            _filtered_spike_times = np.logical_and(
                _spike_times >= time_filter_pre,
                _spike_times < time_filter_post)
            _filtered_spike_rates[period] = \
                np.count_nonzero(_filtered_spike_times) / \
                (current_period_duration * ms)
            for nid in range(all_neurons[pop]):
                _spikes_for_nid = spikes[spikes[:, 0] == nid][:, 1]
                _no_spike_for_nid = np.count_nonzero(np.logical_and(
                    _spikes_for_nid >= time_filter_pre,
                    _spikes_for_nid < time_filter_post))
                per_neuron_spike_count[pop][nid, period] = _no_spike_for_nid
                per_neuron_firing[pop][nid, period] = \
                    _no_spike_for_nid / (current_period_duration * ms)
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
        print("\t{:20}->[{:>8.2f}, {:>8.2f}, {:>8.2f}] Hz".format(
            pop, before, during, after), "per neuron")

    print("=" * 80)
    print("(LaTeX formatting) "
          "Average firing rates before, during and after stimulation "
          "(LaTeX formatting)")
    print("-" * 80)
    for pop in plot_order:
        _x = average_firing_rates[pop] / Hz
        before = _x[0]
        during = _x[1]
        after = _x[2]
        print("\t{:20} & {:>8.2f} & {:>8.2f} & {:>8.2f}".format(
            pop, before, during, after))

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
        _neurons_that_fire = np.sum(
            per_neuron_spike_count[pop], axis=1) > 0
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

        # check that neurons are not both inhibited and excited
        assert np.all(~(np.logical_and(_excited_map, _inhibited_map)))

        excited_filtered_mean = np.mean(_x[_excited_map], axis=0)
        excited_filtered_std = np.std(_x[_excited_map], axis=0)
        excited_before = excited_filtered_mean[0]
        excited_during = excited_filtered_mean[1]
        excited_after = excited_filtered_mean[2]
        excited_before_std = excited_filtered_std[0]
        excited_during_std = excited_filtered_std[1]
        excited_after_std = excited_filtered_std[2]
        print("\t{:20} excited   ->[{:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}] Hz".format(
            pop,
            excited_before, excited_before_std,
            excited_during, excited_during_std,
            excited_after, excited_after_std),
            "per neuron")
        no_excited = np.count_nonzero(_excited_map)
        print("\t\t\t {:6d} excited neurons, i.e. {:7.2%} of cells".format(
            no_excited, no_excited / all_neurons[pop]))

        inhibited_filtered_mean = np.mean(_x[_inhibited_map], axis=0)
        inhibited_filtered_std = np.std(_x[_inhibited_map], axis=0)
        inhibited_before = inhibited_filtered_mean[0]
        inhibited_during = inhibited_filtered_mean[1]
        inhibited_after = inhibited_filtered_mean[2]
        inhibited_before_std = inhibited_filtered_std[0]
        inhibited_during_std = inhibited_filtered_std[1]
        inhibited_after_std = inhibited_filtered_std[2]
        print("\t{:20} inhibited ->[{:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}, {:>8.2f}+-{:>4.1f}] Hz".format(
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
        print("<(LaTeX formatting)>")
        print("{:15} EXCITED  & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f}".format(
            pop,
            excited_before, excited_before_std,
            excited_during, excited_during_std,
            excited_after, excited_after_std))
        print("\t {:6d} ({:7.2%})".format(
            no_excited, no_excited / all_neurons[pop]))
        print("{:15} INHIBITED  & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f}".format(
            pop,
            inhibited_before, inhibited_before_std,
            inhibited_during, inhibited_during_std,
            inhibited_after, inhibited_after_std))
        print("\t {:6d} ({:7.2%})".format(
            no_inhibited, no_inhibited / all_neurons[pop]))
        print("</(LaTeX formatting)>")
        print("-" * 80)
    print("=" * 80)

    if worst_case:
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
        original_conn = np.abs(conn_params[key]["weight"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        # assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95
        _c = Fore.GREEN if is_close else Fore.RED

        print("{:27} -> {}{:4.6f}{} uS".format(
            key, _c, mean, Style.RESET_ALL),
            "c.f. {: 4.6f} uS ({:>7.2%})".format(
                conn_params[key]["weight"], proportion))
    # Report delay values
    print("=" * 80)
    print("Average Delay per projection")
    print("-" * 80)
    for key in final_connectivity:
        conn = conn_dict[key]
        mean = np.mean(conn[:, 3])
        # replace with percentage of difference
        original_conn = np.abs(conn_params[key]["delay"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        # assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95
        _c = Fore.GREEN if is_close else Fore.RED

        print("{:27} -> {}{:4.2f}{} ms".format(
            key, _c, mean, Style.RESET_ALL),
            "c.f. {: 4.2f} ms ({:>7.2%})".format(
                conn_params[key]["delay"], proportion))

    # Report delay values
    print("=" * 80)
    print("(LaTeX formatting) "
          "Average weight per projection "
          "(LaTeX formatting)")
    print("-" * 80)
    for key in final_connectivity:
        conn = conn_dict[key]
        mean = np.mean(conn[:, 2])
        # replace with percentage of difference
        original_conn = np.abs(conn_params[key]["weight"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        # assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95
        _c = Fore.GREEN if is_close else Fore.RED

        print("{:27} & {:4.6f} &".format(key, mean),
              "{: 4.6f} ({:>7.2%})".format(
                  conn_params[key]["weight"], proportion))

    # Check voltage information
    all_voltages = {}
    all_exc_gsyn = {}
    all_inh_gsyn = {}

    if other_recordings is not None:
        print("=" * 80)
        print("Input current analysis")
        print("-" * 80)
        # Looking at the voltage
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
                except AttributeError as ae:
                    all_voltages[pop] = np.squeeze(
                        np.asarray(curr_v.segments[0].analogsignals).T, axis=-1)

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
            nid, tstep = np.unravel_index(np.argmax(v, axis=None), v.shape)
            print("{:20}-> neuron {:>8d} received {:>6d}".format(
                key, int(nid), int(np.max(v))),
                "nA in timestep #{:8d}".format(int(tstep)))
        # Looking at gsyn
        for pop in plot_order:
            try:
                og = other_recordings[pop]
                if og is None:
                    raise KeyError()
                else:
                    curr_exc_gsyn = og['gsyn_exc']
                    curr_inh_gsyn = og['gsyn_inh']
            except KeyError:
                print("No gsyn information for", pop)
                continue

            if (isinstance(curr_exc_gsyn, neo.Block)
                    and isinstance(curr_inh_gsyn, neo.Block)):
                try:
                    all_exc_gsyn[pop] = np.array(curr_exc_gsyn.segments[0].filter(name='gsyn_exc')[0]).T
                    all_inh_gsyn[pop] = np.array(curr_inh_gsyn.segments[0].filter(name='gsyn_inh')[0]).T

                except:
                    all_exc_gsyn[pop] = np.squeeze(
                        np.asarray(curr_exc_gsyn.segments[0].analogsignals).T,
                        axis=-1)
                    all_inh_gsyn[pop] = np.squeeze(
                        np.asarray(curr_inh_gsyn.segments[0].analogsignals).T,
                        axis=-1)

        # report gsyn information
        print("=" * 80)
        print("Excitatory synaptic conductance analysis")
        print("-" * 80)
        for key, g in all_exc_gsyn.items():
            nid, tstep = np.unravel_index(np.argmax(g, axis=None), g.shape)
            print("{:20}-> neuron {:>8d} had an excitatory synaptic "
                  "conductance (g_syn) {:>10.2f}".format(
                key, int(nid), np.max(g)),
                "uS in timestep #{:8d}".format(int(tstep)))
        print("=" * 80)
        print("Inhibitory synaptic conductance analysis")
        print("-" * 80)
        for key, g in all_inh_gsyn.items():
            nid, tstep = np.unravel_index(np.argmax(g, axis=None), g.shape)
            print("{:20}-> neuron {:>8d} had an inhibitory synaptic "
                  "conductance (g_syn) {:>10.2f}".format(
                key, int(nid), np.max(g)),
                "uS in timestep #{:8d}".format(int(tstep)))
    else:
        print("No other recording information present.")

    if worst_case:
        # The following is expensive time wise
        for key, conn in conn_dict.items():
            post_pop = conn_params[key]["post"]
            pre_pop = conn_params[key]["pre"]
            curr_spikes = all_spikes[pre_pop]
            for nid, t in curr_spikes:
                nid = int(nid)
                times = int(t * time_to_bin_conversion)
                targets = conn[conn[:, 0] == nid][:, 1].astype(int)
                inc_spike_count[post_pop][targets, times] += 1

    if conn_exists and worst_case:
        print("=" * 80)
        print("Incoming spikes statistics")
        print("-" * 80)
        for pop, counts in inc_spike_count.items():
            nid, tstep = np.unravel_index(np.argmax(counts, axis=None), counts.shape)
            print("{:20}-> neuron {:>8d} saw {:>6d}".format(
                pop, int(nid), int(np.max(counts))),
                "spikes in timestep #{:8d}".format(int(tstep)))
            # Print worst case statistic for the population
            maxs = np.max(counts, axis=1)
            assert maxs.size == all_neurons[pop]
            print("\t# spikes: mean {:8.4f}, "
                  "mode {:12.8f}, "
                  "std {:8.4f}".format(
                np.mean(maxs), scipy.stats.mode(maxs).mode[0], np.std(maxs)
            ))

    print("=" * 80)
    print("Plotting figures...")
    print("-" * 80)

    wanted_times = np.linspace(0, (simtime / ms), 6).astype(int)

    common_highlight_values = {
        'start': stim_period_start,
        'stop': stim_period_end,
        'increment': timestep / ms,
    }

    common_values_for_plots = {
        'plot_order': plot_order,
        'wanted_times': wanted_times,
        'time_to_bin_conversion': time_to_bin_conversion,
        'fig_folder': sim_fig_folder,
        'highlight_stim': highlight_stim,
        'common_highlight_values': common_highlight_values,
    }

    plot_analog_signal(all_exc_gsyn, variable_name="gsyn_exc",
                       ylabel="Exc synaptic conductance ($\mu S$)",
                       **common_values_for_plots)
    plot_analog_signal(all_inh_gsyn, variable_name="gsyn_inh",
                       ylabel="Inh synaptic conductance ($\mu S$)",
                       **common_values_for_plots)
    plot_analog_signal(all_voltages, variable_name="v",
                       ylabel="Membrane potential (mV)",
                       **common_values_for_plots)

    # Plot distribution of worst case spikes per population
    if conn_exists and worst_case:
        print("Plotting histogram of worst spike counts")
        for index, pop in enumerate(plot_order):
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            counts = inc_spike_count[pop]
            maxs = np.max(counts, axis=1)
            assert maxs.size == all_neurons[pop]
            plt.hist(maxs, bins=20, color=viridis_cmap(index / (n_plots + 1)),
                     rasterized=True,
                     edgecolor='k')

            plt.title(use_display_name(pop))

            plt.ylabel("Count")
            plt.xlabel("Max # of spikes per neuron")
            plt.tight_layout()
            save_figure(
                plt,
                os.path.join(sim_fig_folder,
                             "max_spikes_per_neuron_in_{}").format(pop),
                extensions=[".png", ".pdf"])
            plt.close(f)

    # raster plot including ALL populations
    print("Plotting spiking raster plot for all populations")
    f, axes = plt.subplots(len(all_spikes.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=400)
    for index, pop in enumerate(plot_order):
        curr_ax = axes[index]
        # spike raster
        _times = all_spikes[pop][:, 1]
        _ids = all_spikes[pop][:, 0]
        if highlight_stim:
            highlight_area(curr_ax, pop, **common_highlight_values)

        curr_ax.scatter(_times,
                        _ids,
                        color=viridis_cmap(index / (n_plots + 1)),
                        s=.5, rasterized=True)
        curr_ax.set_title(use_display_name(pop))
    plt.xlabel("Time (ms)")
    # plt.suptitle((use_display_name(simulator)+"\n")
    f.tight_layout()
    save_figure(plt, os.path.join(sim_fig_folder, "raster_plots"),
                extensions=['.png', '.pdf'])
    plt.close(f)

    # raster plot + PSTH including ALL populations
    print("Plotting spiking raster plot + PSTH for all populations")
    f, axes = plt.subplots(2 * len(all_spikes.keys()), 1,
                           figsize=(14, 30), sharex=True, dpi=500)
    for index, pop in enumerate(np.repeat(plot_order, 2)):
        curr_ax = axes[index]
        if highlight_stim:
            highlight_area(curr_ax, pop, **common_highlight_values)
        if index % 2 == 0:
            # spike raster
            _times = all_spikes[pop][:, 1]
            _ids = all_spikes[pop][:, 0]
            curr_ax.scatter(_times,
                            _ids,
                            color=viridis_cmap(int(index / 2) / (n_plots + 1)),
                            s=.5, rasterized=True)
            curr_ax.set_title(use_display_name(pop))
            curr_ax.set_ylabel("NID")
        else:
            curr_ax.bar(np.arange(spikes_per_timestep[pop].size) * timestep / ms,
                        spikes_per_timestep[pop],
                        color=viridis_cmap(int(index / 2) / (n_plots + 1)),
                        rasterized=True)
            curr_ax.set_ylabel("Count")

    plt.xlabel("Time (ms)")
    # plt.suptitle((use_display_name(simulator) + "\n")
    f.tight_layout()
    save_figure(plt, os.path.join(sim_fig_folder, "raster_and_psth_plots"),
                extensions=['.png', ])
    plt.close(f)

    # raster + PSTH for each population
    print("Plotting spiking raster plot + PSTH for each population")
    for index, pop in enumerate(plot_order):
        print("\t{:20}".format(pop), end=' ')
        f, (ax_0, ax_1) = plt.subplots(2, 1, figsize=(9, 9),
                                       sharex=True, dpi=400)

        # spike raster
        _times = all_spikes[pop][:, 1]
        _ids = all_spikes[pop][:, 0]

        if highlight_stim:
            highlight_area(ax_0, pop, **common_highlight_values)
        ax_0.scatter(_times,
                     _ids,
                     color=viridis_cmap(index / (n_plots + 1)),
                     s=.5, rasterized=True)
        ax_0.set_ylabel("NID")

        # PSTH
        if highlight_stim:
            highlight_area(ax_1, pop, **common_highlight_values)
        ax_1.bar(np.arange(spikes_per_timestep[pop].size) * timestep / ms,
                 spikes_per_timestep[pop],
                 color=viridis_cmap(index / (n_plots + 1)))
        ax_1.set_ylabel("Count")
        plt.xlabel("Time (ms)")
        save_figure(plt, os.path.join(sim_fig_folder,
                                      "{}_raster_and_psth".format(pop)),
                    extensions=['.png', ])
        plt.close(f)
        print("SUCCESS")

    # raster + PSTH + voltage for each population
    print("Plotting spiking raster plot + PSTH + voltage for each population")
    for index, pop in enumerate(plot_order):
        # plot voltage traces
        print("\t{:20}".format(pop), end=' ')
        if pop in ["glomerulus", "mossy_fibers"]:
            print("FAIL -- spike source")
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            plt.close(f)
            continue
        f, (ax_0, ax_1, ax_2) = plt.subplots(3, 1, figsize=(9, 12),
                                             sharex=True, dpi=400)

        try:
            # spike raster
            if highlight_stim:
                highlight_area(ax_0, pop, **common_highlight_values)
            ax_0.scatter(all_spikes[pop][:, 1],
                         all_spikes[pop][:, 0],
                         color=viridis_cmap(index / (n_plots + 1)),
                         s=.5, rasterized=True)
            ax_0.set_ylabel("NID")
            # PSTH
            if highlight_stim:
                highlight_area(ax_1, pop, **common_highlight_values)
            ax_1.bar(np.arange(spikes_per_timestep[pop].size) * timestep / ms,
                     spikes_per_timestep[pop],
                     color=viridis_cmap(index / (n_plots + 1)),
                     rasterized=True)
            ax_1.set_ylabel("Count")

            # voltage
            pop_exc_g = all_voltages[pop]
            for _ind, _trace in enumerate(pop_exc_g):
                ax_2.plot(np.arange(_trace.size) * timestep / ms, _trace,
                          color=viridis_cmap(index / (n_plots + 1)))

            ax_2.set_ylabel("Membrane potential (mV)")
            plt.xlabel("Time (ms)")

            save_figure(plt, os.path.join(sim_fig_folder,
                                          "{}_raster_psth_and_voltage".format(pop)),
                        extensions=['.png', ])
            plt.close(f)
            print("SUCCESS")
        except:
            print("FAIL")
    # plot .1 ms PSTH
    print("Plotting PSTH for each timestep")
    f, axes = plt.subplots(len(spikes_per_timestep.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=400)
    for index, pop in enumerate(plot_order):
        if highlight_stim:
            highlight_area(axes[index], pop, **common_highlight_values)
        axes[index].bar(np.arange(spikes_per_timestep[pop].size),
                        spikes_per_timestep[pop],
                        color=viridis_cmap(index / (n_plots + 1)),
                        rasterized=True)
        axes[index].set_title(use_display_name(pop))
    plt.xticks(wanted_times * time_to_bin_conversion, wanted_times)
    # plt.suptitle((use_display_name(simulator) + "\n")

    save_figure(plt, os.path.join(sim_fig_folder,
                                  "timestep_psth"),
                extensions=['.png', ])
    plt.close(f)

    # plot 3 ms PSTH
    print("Plotting PSTH in bins of 3 ms")
    f, axes = plt.subplots(len(spikes_per_3ms.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=400)
    for index, pop in enumerate(plot_order):
        axes[index].bar(np.arange(spikes_per_3ms[pop].size), spikes_per_3ms[pop],
                        color=viridis_cmap(index / (n_plots + 1)),
                        rasterized=True)
        axes[index].set_title(use_display_name(pop))

    plt.xticks(wanted_times * time_to_bin_conversion / bins_in_3ms, wanted_times)
    # plt.suptitle((use_display_name(simulator) + "\n")

    save_figure(plt, os.path.join(sim_fig_folder,
                                  "timestep_psth_3ms"),
                extensions=['.png', '.pdf'])
    plt.close(f)

    # plot firing rate histogram per PSTH region
    print("Plotting firing rate histograms")
    f, axes = plt.subplots(len(plot_order), 3,
                           figsize=(14, 20), sharex=True, dpi=400)
    for index, pop in enumerate(plot_order):
        for period in range(stimulus_periods):
            curr_ax = axes[index, period]
            curr_ax.hist(per_neuron_firing[pop][:, period],
                         color=viridis_cmap(index / (n_plots + 1)),
                         bins=20, rasterized=True)
            if period == 1:
                curr_ax.set_title(use_display_name(pop))
            curr_ax.set_xlabel("Hz")
            curr_ax.xaxis.set_tick_params(which='both',
                                          labelbottom=True)
            curr_ax.set_xticks([0, 75, 150])

    f.tight_layout()

    save_figure(plt, os.path.join(sim_fig_folder,
                                  "neuron_firing_rate_hist"),
                extensions=['.png', '.pdf'])
    plt.close(f)
    # TODO plot weight histogram

    # TODO plot centred connectivity
    print("=" * 80)


def compare_results(file_1, file_2, fig_folder, dark_background):
    if dark_background:
        plt.style.use('dark_background')

    # Retrieve results file
    try:
        data_1 = np.load(file_1, allow_pickle=True)
    except FileNotFoundError:
        data_1 = np.load(file_1 + ".npz", allow_pickle=True)
        file_1 += ".npz"

    # Retrieve results file
    try:
        data_2 = np.load(file_2, allow_pickle=True)
    except FileNotFoundError:
        data_2 = np.load(file_2 + ".npz", allow_pickle=True)
        file_2 += ".npz"

    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)

    # Create figures folder for this results_file
    base_name_file_1 = str(ntpath.basename(file_1))[:-4]
    base_name_file_2 = str(ntpath.basename(file_2))[:-4]

    sim_fig_folder = os.path.join(fig_folder,
                                  base_name_file_1 + "_vs_" + base_name_file_2)
    if not os.path.isdir(sim_fig_folder) and not os.path.exists(sim_fig_folder):
        os.mkdir(sim_fig_folder)

    # Set up colours
    color_init(strip=False)

    # Plotting results for ...
    print("=" * 80)
    print("Plotting comaprison results between", file_1, "and", file_2)
    print("-" * 80)

    all_spikes_1 = data_1['all_spikes'].ravel()[0]
    all_spikes_2 = data_2['all_spikes'].ravel()[0]

    other_recordings_1 = data_1['other_recordings'].ravel()[0]
    other_recordings_2 = data_2['other_recordings'].ravel()[0]

    # Compute plot order
    plot_order = get_plot_order(all_spikes_1.keys())
    n_plots = float(len(plot_order))

    sim_params_1 = data_1['simulation_parameters'].ravel()[0]
    sim_params_2 = data_2['simulation_parameters'].ravel()[0]
    timestep = sim_params_1['argparser']['timestep'] * ms

    # TODO assert that parameters are the same between the simulations

    simtime = data_1['simtime'] * pq.ms
    assert (simtime == (data_2['simtime'] * pq.ms))

    simulators = [sim_params_1['argparser']['simulator'],
                  sim_params_2['argparser']['simulator']]

    reporting_format_string = "{:30}-{:30}:{:6.3f}+-{:6.3f}"

    f = plt.figure(1, figsize=(9, 9), dpi=400)
    plt.close(f)

    all_corr_coef = {}
    all_distances = {}
    all_tiling_coef = {}
    all_isi = {}
    all_cv = {}
    all_cv2 = {}
    print("=" * 80)
    print("Spike train analysis (avg +- std)")
    # https://elephant.readthedocs.io/en/latest/modules.html
    print("-" * 80)
    for pop in plot_order:
        pop_1_spikes = all_spikes_1[pop].segments[0].spiketrains
        pop_2_spikes = all_spikes_2[pop].segments[0].spiketrains

        corr_coef_per_neuron = []
        spike_time_tiling_coef_per_neuron = []
        van_rossum_dists = []
        isis = {0: [],
                1: []}
        cvs = {0: [],
               1: []}
        cv2s = {0: [],
                1: []}
        # I guess we need to look at each neuron?
        for (p1, p2) in zip(pop_1_spikes, pop_2_spikes):
            cc_matrix = elephant.spike_train_correlation.corrcoef(
                elephant.conversion.BinnedSpikeTrain([p1, p2], binsize=5 * pq.ms),
                fast=True)
            corr_coef_per_neuron.append(cc_matrix[0, 1])

            tiling = elephant.spike_train_correlation.spike_time_tiling_coefficient(
                p1, p2, dt=1 * pq.ms)
            spike_time_tiling_coef_per_neuron.append(tiling)

            dist = elephant.spike_train_dissimilarity.van_rossum_dist(
                [p1, p2], tau=1 * pq.ms)[0, 1]
            van_rossum_dists.append(dist)

            # Inter spike intervals
            isis[0].append(np.array(elephant.statistics.isi(p1 / pq.ms)))
            isis[1].append(np.array(elephant.statistics.isi(p2 / pq.ms)))

            # Coefficient of variation
            cvs[0].append(elephant.statistics.cv(p1))
            cvs[1].append(elephant.statistics.cv(p2))

            # CV2 https://elephant.readthedocs.io/en/latest/reference/toctree/statistics/elephant.statistics.cv2.html#elephant.statistics.cv2
            cv2s[0].append(elephant.statistics.cv2(p1, with_nan=True))
            cv2s[1].append(elephant.statistics.cv2(p2, with_nan=True))

        print(reporting_format_string.format(
            pop, "corr coef",
            np.nanmean(corr_coef_per_neuron),
            np.nanstd(corr_coef_per_neuron)))
        print(reporting_format_string.format(
            "", "tiling coef",
            np.nanmean(spike_time_tiling_coef_per_neuron),
            np.nanstd(spike_time_tiling_coef_per_neuron)))
        print(reporting_format_string.format(
            "", "van rossum dist",
            np.nanmean(van_rossum_dists),
            np.nanstd(van_rossum_dists)))

        all_corr_coef[pop] = corr_coef_per_neuron
        all_distances[pop] = van_rossum_dists
        all_tiling_coef[pop] = spike_time_tiling_coef_per_neuron

        all_isi[pop] = isis
        all_cv[pop] = cvs
        all_cv2[pop] = cv2s

    for k, v in all_isi.items():
        for dk, dv in v.items():
            np_of_dv = np.array(dv)
            try:
                flatten_list = np.hstack(np_of_dv).ravel()
                all_isi[k][dk] = flatten_list
            except:
                pass

    all_coherence = {}
    all_freqs = {}
    all_lags = {}
    if other_recordings_1 is not None and other_recordings_2 is not None:
        print("=" * 80)
        print("Analogue value analysis (avg +- std)")
        # https://elephant.readthedocs.io/en/latest/reference/spectral.html
        print("-" * 80)
        for signal in ['v', 'gsyn_exc', 'gsyn_inh']:
            all_coherence[signal] = {}
            all_freqs[signal] = {}
            all_lags[signal] = {}
            for pop in plot_order:
                if pop not in other_recordings_1.keys():
                    continue
                pop_1_v = other_recordings_1[pop][signal].segments[0].analogsignals
                pop_2_v = other_recordings_2[pop][signal].segments[0].analogsignals

                welch_cohere_freqs_per_neuron = []
                welch_cohere_lags_per_neuron = []
                welch_cohere_coher_per_neuron = []
                # I guess we need to look at each neuron?
                for (p1, p2) in zip(pop_1_v, pop_2_v):
                    if len(p1) != len(p2):
                        min_len = min(len(p1), len(p2))
                        p1 = p1[:min_len]
                        p2 = p2[:min_len]
                    f, c, l = elephant.spectral.welch_cohere(
                        p1, p2)
                    welch_cohere_freqs_per_neuron.append(f)
                    welch_cohere_lags_per_neuron.append(l)
                    welch_cohere_coher_per_neuron.append(c)

                welch_cohere_freqs_per_neuron = np.asarray(welch_cohere_freqs_per_neuron)
                welch_cohere_lags_per_neuron = np.asarray(welch_cohere_lags_per_neuron)
                welch_cohere_coher_per_neuron = np.asarray(welch_cohere_coher_per_neuron)
                print(reporting_format_string.format(
                    pop + " " + signal, "coherence ([0, 1])",
                    np.nanmean(welch_cohere_coher_per_neuron),
                    np.nanstd(welch_cohere_coher_per_neuron))
                )
                print("{:30}-{:30}:{:6.3f}+-{:6.3f}".format(
                    "", "lag (phase, [-pi, +pi])",
                    np.nanmean(welch_cohere_lags_per_neuron),
                    np.nanstd(welch_cohere_lags_per_neuron)))
                all_coherence[signal][pop] = welch_cohere_coher_per_neuron
                all_lags[signal][pop] = welch_cohere_coher_per_neuron
                all_freqs[signal][pop] = welch_cohere_freqs_per_neuron

    print("=" * 80)
    print("Plotting figures...")
    print("-" * 80)

    # raster plot including ALL populations
    print("Plotting spiking raster plot for all populations")

    wanted_times = np.linspace(0, (simtime / ms), 6).astype(int)

    common_highlight_values = {
        'start': None,  # stim_period_start,
        'stop': None,  # stim_period_end,
        'increment': timestep / ms,
    }

    common_values_for_plots = {
        'plot_order': plot_order,
        'wanted_times': None,  # wanted_times,
        'time_to_bin_conversion': None,  # time_to_bin_conversion,
        'fig_folder': sim_fig_folder,
        'highlight_stim': None,  # highlight_stim,
        'common_highlight_values': None,  # common_highlight_values,
    }

    # Plotting boxplots
    plot_boxplot_for_all_pops(all_corr_coef, "corr_coeff",
                              xlabel="Population",
                              ylabel="Correlation Coefficient",
                              **common_values_for_plots)

    plot_boxplot_for_all_pops(all_distances, "van_rostum_dist",
                              xlabel="Population",
                              ylabel="van Rostum Distance",
                              **common_values_for_plots)

    plot_boxplot_for_all_pops(all_tiling_coef, "tiling_coeff",
                              xlabel="Population",
                              ylabel="Tiling Coefficient",
                              **common_values_for_plots)
    # Plot side by side boxplots
    plot_sbs_boxplot_for_all_pops(all_isi, "isi",
                                  xlabel="Population",
                                  ylabel="ISI (ms)",titles=simulators,
                                  **common_values_for_plots)

    plot_sbs_boxplot_for_all_pops(all_cv, "cv",
                                  xlabel="Population",
                                  ylabel="CV",titles=simulators,
                                  **common_values_for_plots)

    plot_sbs_boxplot_for_all_pops(all_cv2, "cv2",
                                  xlabel="Population",
                                  ylabel="CV2",titles=simulators,
                                  **common_values_for_plots)

    # Plot side by side histograms
    plot_sbs_histogram_for_all_pops(all_isi, "isi",
                                    xlabel="ISI",
                                    ylabel="count", titles=simulators,
                                    **common_values_for_plots)

    plot_sbs_histogram_for_all_pops(all_cv, "cv",
                                    xlabel="CV",
                                    ylabel="count", titles=simulators,
                                    **common_values_for_plots)

    plot_sbs_histogram_for_all_pops(all_cv2, "cv2",
                                    xlabel="CV2",
                                    ylabel="count", titles=simulators,
                                    **common_values_for_plots)

    # Plotting histograms
    plot_histogram_for_all_pops(all_corr_coef, "corr_coeff",
                                xlabel="Correlation Coefficient",
                                ylabel="count",
                                **common_values_for_plots)

    plot_histogram_for_all_pops(all_distances, "van_rostum_dist",
                                xlabel="van Rostum Distance",
                                ylabel="count",
                                **common_values_for_plots)

    plot_histogram_for_all_pops(all_tiling_coef, "tiling_coeff",
                                xlabel="Tiling Coefficient",
                                ylabel="count",
                                **common_values_for_plots)

    print("=" * 80)


def plot_sbs_histogram_for_all_pops(data, variable_name,
                                    xlabel, ylabel, plot_order,
                                    wanted_times, time_to_bin_conversion,
                                    fig_folder,
                                    highlight_stim, common_highlight_values,
                                    titles):
    print("Plotting SIDE BY SIDE histogram for {} for all populations".format(
        variable_name))
    f, axes = plt.subplots(len(data.keys()), 2,
                           figsize=(14, 20), sharex='row', sharey='row', dpi=400)
    n_plots = len(plot_order)
    minimus = 0
    maximus = 0
    for index, pop in enumerate(plot_order):
        curr_ax = axes[index][0]
        curr_ax2 = axes[index][1]
        curr_pop = data[pop]
        curr_ax.hist(np.asarray(curr_pop[0]).ravel(), bins=21, edgecolor='k',
                     color=color_for_index(index, n_plots),
                     rasterized=True, normed=True)
        curr_ax2.hist(np.asarray(curr_pop[1]).ravel(), bins=21, edgecolor='k',
                      color=color_for_index(index, n_plots),
                      rasterized=True, normed=True)
        if index == 0:
            curr_ax.set_title(use_display_name(titles[0]))
            curr_ax2.set_title(use_display_name(titles[1]))
        curr_ax.set_ylabel(use_display_name(pop) + "\n" + ylabel)
        curr_ax.grid(True, which="major", axis="both")
        curr_ax2.grid(True, which="major", axis="both")

        minimus = min(minimus, min(data[pop][0]), min(data[pop][1]))
        maximus = max(maximus, max(data[pop][0]), max(data[pop][1]))

    plt.xlabel(xlabel)
    plt.xlim([minimus - (.1 * maximus), 1.1 * maximus])
    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "sbs_all_pop_{}".format(variable_name)),
                extensions=['.png', '.pdf'])
    plt.close(f)


def plot_histogram_for_all_pops(data, variable_name, xlabel, ylabel, plot_order,
                                wanted_times, time_to_bin_conversion,
                                fig_folder,
                                highlight_stim, common_highlight_values):
    print("Plotting histogram for {} for all populations".format(variable_name))
    f, axes = plt.subplots(len(data.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=400)
    n_plots = len(plot_order)
    minimus = 0
    maximus = 0
    for index, pop in enumerate(plot_order):
        curr_ax = axes[index]
        curr_ax.hist(data[pop], bins=21, edgecolor='k',
                     color=color_for_index(index, n_plots),
                     rasterized=True, normed=True)
        curr_ax.set_title(use_display_name(pop))
        curr_ax.set_ylabel(ylabel)
        curr_ax.grid(True, which="major", axis="both")

        minimus = min(minimus, min(data[pop]))
        maximus = max(maximus, max(data[pop]))
    plt.xlabel(xlabel)
    plt.xlim([minimus - (.1 * maximus), 1.1 * maximus])
    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "all_pop_{}".format(variable_name)),
                extensions=['.png', '.pdf'])
    plt.close(f)


def plot_boxplot_for_all_pops(data, variable_name, xlabel, ylabel, plot_order,
                              wanted_times, time_to_bin_conversion,
                              fig_folder,
                              highlight_stim, common_highlight_values):
    print("Plotting boxplot for {} for all populations".format(variable_name))
    assembled_results = []
    for index, pop in enumerate(plot_order):
        assembled_results.append(data[pop])
    n_plots = len(plot_order)
    # stlye the median of boxplots
    medianprops = dict(color='#414C82', linewidth=1.5)

    f = plt.figure(figsize=(12, 8), dpi=600)
    for index, pop in enumerate(plot_order):
        plt.boxplot(data[pop], notch=True, positions=[index + 1],
                    medianprops=dict(color=color_for_index(index, n_plots), linewidth=1.5))
    # plt.boxplot(assembled_results, notch=True, medianprops=medianprops)
    # plt.title(use_display_name(pop))
    plt.ylabel(ylabel)
    plt.xlim([0, n_plots + 1])
    plt.grid(True, which="major", axis="y")
    plt.xlabel(xlabel)
    xtick_display_names = [use_display_name(x) for x in plot_order]
    _, labels = plt.xticks(np.arange(len(plot_order)) + 1, xtick_display_names)

    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "bp_all_pop_{}".format(variable_name)),
                extensions=['.png', '.pdf'])
    plt.close(f)


def plot_sbs_boxplot_for_all_pops(data, variable_name, xlabel, ylabel,
                                  plot_order,
                                  wanted_times, time_to_bin_conversion,
                                  fig_folder,
                                  highlight_stim, common_highlight_values,
                                  titles):
    print("Plotting SIDE BY SIDE boxplot for {} for all populations".format(
        variable_name))
    assembled_results = []
    for index, pop in enumerate(plot_order):
        assembled_results.append(data[pop])
    n_plots = len(plot_order)
    spacer = 4
    bp_width = 0.7
    f = plt.figure(figsize=(12, 8), dpi=600)
    for index, pop in enumerate(plot_order):
        curr_pop = data[pop]
        bp0 = plt.boxplot(np.asarray(curr_pop[0]).ravel(), notch=True,
                    positions=[spacer * index + 1],
                    medianprops=dict(color=viridis_cmap(0.0),
                                     linewidth=1.5),
                    widths=bp_width)
        bp1 = plt.boxplot(np.asarray(curr_pop[1]).ravel(), notch=True,
                    positions=[spacer * index + 2],
                    medianprops=dict(color=viridis_cmap(.9),
                                     linewidth=1.5),
                    widths=bp_width)
    title_display_names = [use_display_name(x) for x in titles]
    plt.legend([bp0["medians"][0], bp1["medians"][0]], title_display_names,
               loc='best')
    plt.ylabel(ylabel)
    plt.xlim([0, spacer * n_plots])
    plt.grid(True, which="major", axis="y")
    plt.xlabel(xlabel)
    xtick_display_names = [use_display_name(x) for x in plot_order]
    _, labels = plt.xticks(spacer * np.arange(len(plot_order)) + 1.5,
                           xtick_display_names)

    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "bp_sbs_all_pop_{}".format(variable_name)),
                extensions=['.png', '.pdf'])
    plt.close(f)


if __name__ == "__main__":
    from spinncer.analysis_argparser import *

    if analysis_args.compare and len(analysis_args.compare) > 0:
        if len(analysis_args.compare) % 2 == 1:
            raise ValueError("The number of archives to compare is meant to "
                             "be a multiple of 2.")
        for i in range(0, len(analysis_args.compare), 2):
            compare_results(analysis_args.compare[i],
                            analysis_args.compare[i + 1],
                            analysis_args.figures_dir,
                            dark_background=analysis_args.dark_background)
    if analysis_args.input and len(analysis_args.input) > 0:
        for in_file in analysis_args.input:
            spike_analysis(in_file, analysis_args.figures_dir,
                           worst_case=analysis_args.worst_case_spikes,
                           delay_sensitive=analysis_args.consider_delays,
                           dark_background=analysis_args.dark_background,
                           highlight_stim=analysis_args.highlight_stim)
