from spinncer.analysis_common import *
import traceback

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
                       highlight_stim, common_highlight_values,
                       scale_xticks=1):
    print("Plotting {} traces for each population".format(variable_name))
    for index, pop in enumerate(plot_order):
        try:
            values_for_pop = data[pop]
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            for _ind, _trace in enumerate(values_for_pop):
                plt.plot(_trace,
                         color=color_for_index(_ind, values_for_pop.shape[0]),
                         rasterized=True)
            plt.xticks(wanted_times * time_to_bin_conversion * scale_xticks,
                       wanted_times)
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


def unpack_abcd(sh_exc, sh_exc_2):
    a = ((sh_exc) & 0xFFFF)
    b = ((sh_exc >> 16) & 0xFFFF)
    c = ((sh_exc_2) & 0xFFFF)
    d = ((sh_exc_2 >> 16) & 0xFFFF)
    return a, b, c, d


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
    for cell_name in plot_order:
        params = cell_params[cell_name]
        if cell_name == "glomerulus":
            continue
        print("{} & {} & {} & {} & {} & {} & {} & {} & {} & {} \\\\".format(
            use_display_name(cell_name),
            params["cm"],
            params["i_offset"],
            params["tau_m"],
            params["tau_refrac"],
            params["tau_syn_E"],
            params["tau_syn_I"],
            params["v_reset"],
            params["v_rest"],
            params["v_thresh"],
        )
        )
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
    # Number of post-synaptic hits in each timestep
    all_post_hits = {}
    all_max_exc_spikes = {}
    all_max_inh_spikes = {}
    all_exc_spike_counts = {}
    all_inh_spikes_counts = {}
    print("=" * 80)
    print("Maximum number of generated spikes per timestep")
    print("-" * 80)

    # Elephant computed metrics
    elephant_psths = {}
    elephant_instantaneous_rates = {}
    elephant_timestep = sim_params['argparser']['timestep'] * pq.ms
    elephant_simtime = data['simtime'] * pq.ms
    try:
        for pop in plot_order:
            curr_spikes = neo_all_spikes[pop].segments[0].spiketrains
            curr_inst_rates = \
                elephant.statistics.instantaneous_rate(
                    curr_spikes,
                    sampling_period=elephant_timestep,
                    t_start=0 * pq.ms,
                    t_stop=elephant_simtime
                )
            curr_psths = \
                elephant.statistics.time_histogram(
                    curr_spikes,
                    binsize=elephant_timestep,
                    t_start=0 * pq.ms,
                    t_stop=elephant_simtime
                )
            elephant_psths[pop] = curr_psths.ravel()
            # Save the average instantaneous rate
            elephant_instantaneous_rates[pop] = curr_inst_rates.ravel() / all_neurons[pop]
    except:
        pass

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
             (spikes_per_timestep[pop].size - no_timesteps)),
            # This should be 0 or 1 corresponding to SpiNNaker or NEST
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
        print("{:15} EXCITED {:6d} ({:7.2%}) & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f}".format(
            pop,
            int(no_excited), no_excited / all_neurons[pop],
            excited_before, excited_before_std,
            excited_during, excited_during_std,
            excited_after, excited_after_std))

        print("{:15} INHIBITED {:6d} ({:7.2%}) & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f} & "
              "{:>8.2f}$\pm${:>4.1f}".format(
            pop,
            int(no_inhibited), no_inhibited / all_neurons[pop],
            inhibited_before, inhibited_before_std,
            inhibited_during, inhibited_during_std,
            inhibited_after, inhibited_after_std))
        print("</(LaTeX formatting)>")
        print("-" * 80)
    print("=" * 80)

    if worst_case:
        # Count incoming spikes only if we care -- this takes a while
        inc_spike_count = {k: np.zeros((all_neurons[k], no_timesteps + 1))
                           for k in all_neurons.keys()}

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
    print("Connection Name ")
    print("{:27} | {:10} | ".format("Connection Name", "Def. W"),
          "{:20}".format("SpiNN. W"))
    for key in final_connectivity:
        conn = conn_dict[key]
        mean = np.mean(conn[:, 2])
        # replace with percentage of difference
        original_conn = np.abs(conn_params[key]["weight"])
        orig_conn_sign = np.sign(conn_params[key]["weight"])
        if mean < original_conn:
            proportion = mean / original_conn
        else:
            proportion = original_conn / mean
        diff = original_conn - mean
        prop_diff = orig_conn_sign * (diff / original_conn)
        # assert (0 <= proportion <= 1), proportion
        is_close = proportion >= .95

        print("{:27} & {:4.10f} &".format(key, conn_params[key]["weight"]),
              "{: 4.10f} ({:>7.2%})".format(
                  mean, prop_diff))

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
        # Report statistics here
        for key, v in all_voltages.items():
            nid, tstep = np.unravel_index(np.argmax(v, axis=None), v.shape)
            print("{:20}-> neuron {:>8d} received {:>6d}".format(
                key, int(nid), int(np.max(v))),
                "nA in timestep #{:8d}".format(int(tstep)))
            # THIS IS BROKEN! it will be removed soon
            # # Also treat voltage as if it's a piggybacked value packaging
            # # counts of number of post-synaptic hits for spikes
            # max_v = np.max(v, axis=0)
            # max_g = np.max(all_exc_gsyn[key], axis=0)
            # abcd_view = max_v.astype(np.uint32).view(dtype=[
            #     #                                   ('d', np.uint16),
            #     #                                   ('c', np.uint16),
            #     ('b', np.uint16),
            #     ('a', np.uint16)
            # ])
            # abcd_view_part_2 = max_g.astype(np.uint32).view(dtype=[
            #     ('d', np.uint16),
            #     ('c', np.uint16),
            #     #                                   ('b', np.uint16),
            #     #                                   ('a', np.uint16)
            # ])
            # pd_view = pd.DataFrame(abcd_view)
            # pd_view_part_2 = pd.DataFrame(abcd_view_part_2)
            # pd_view = pd.concat([pd_view, pd_view_part_2], axis=1, sort=False)
            #
            # pd_view.describe()
            # # this is meaningless using default tools
            # all_post_hits[key] = pd_view
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
            # this is meaningless using default tools
            all_max_exc_spikes[key] = np.max(g, axis=0)
        print("=" * 80)
        print("Inhibitory synaptic conductance analysis")
        print("-" * 80)
        for key, g in all_inh_gsyn.items():
            nid, tstep = np.unravel_index(np.argmax(g, axis=None), g.shape)
            print("{:20}-> neuron {:>8d} had an inhibitory synaptic "
                  "conductance (g_syn) {:>10.2f}".format(
                key, int(nid), np.max(g)),
                "uS in timestep #{:8d}".format(int(tstep)))
            # this is meaningless using default tools
            all_max_inh_spikes[key] = np.max(g, axis=0)
    else:
        print("No other recording information present.")

    if worst_case:
        per_conn_worst_spikes = {}
        for k in all_neurons.keys():
            per_conn_worst_spikes[k] = {}
        # The following is expensive time wise
        for key, conn in conn_dict.items():
            post_pop = conn_params[key]["post"]
            pre_pop = conn_params[key]["pre"]
            curr_spikes = all_spikes[pre_pop]
            inc_spikes_for_this_conn = np.zeros(inc_spike_count[post_pop].shape)
            for nid, t in curr_spikes:
                nid = int(nid)
                times = int(t * time_to_bin_conversion)
                targets = conn[conn[:, 0] == nid][:, 1].astype(int)
                inc_spikes_for_this_conn[targets, times] += 1
            inc_spike_count[post_pop] += inc_spikes_for_this_conn
            per_conn_worst_spikes[post_pop][key] = inc_spikes_for_this_conn

    if conn_exists and worst_case:
        print("=" * 80)
        print("Incoming spikes statistics")
        print("-" * 80)
        for pop in plot_order:
            counts = inc_spike_count[pop]
            nid, tstep = np.unravel_index(np.argmax(counts, axis=None), counts.shape)
            print("{:20}-> neuron {:>8d} saw {:>6d}".format(
                pop, int(nid), int(np.max(counts))),
                "spikes in timestep #{:8d}".format(int(tstep)))
            # Print worst case statistic for the population
            maxs = np.max(counts, axis=1)
            assert maxs.size == all_neurons[pop]
            print("\t# spikes: mean {:8.4f}, "
                  "std {:8.4f}".format(
                np.mean(maxs), np.std(maxs)
            ))
            for conn_key in per_conn_worst_spikes[pop]:
                curr_nid_inc_spikes = per_conn_worst_spikes[pop][conn_key][nid]
                print("\t\t{:10} contributed {:8d} spikes".format(
                    conn_key, int(curr_nid_inc_spikes[tstep])
                ))

    print("=" * 80)
    print("Plotting figures...")
    print("-" * 80)

    wanted_times = np.linspace(0, (simtime / ms), 6).astype(int)
    stim_wanted_times = np.linspace(time_filter[1] - 50,
                                    time_filter[2] + 50, 4).astype(int)
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

    # plot .1 ms PSTH
    print("Plotting elephant Instantaneous rates")
    try:
        f, axes = plt.subplots(len(spikes_per_timestep.keys()), 1,
                               figsize=(14, 20), sharex=True, dpi=400)
        for index, pop in enumerate(plot_order):
            if highlight_stim:
                highlight_area(axes[index], pop, **common_highlight_values)
            curr_psth = elephant_instantaneous_rates[pop]
            axes[index].bar(np.arange(curr_psth.size),
                            curr_psth,
                            color=viridis_cmap(index / (n_plots + 1)),
                            rasterized=True)
            axes[index].set_title(use_display_name(pop))
            axes[index].set_ylabel("Hz")
        plt.xticks(wanted_times * time_to_bin_conversion, wanted_times)
        plt.xlabel("Time (ms)")
        save_figure(plt, os.path.join(sim_fig_folder,
                                      "elephant_instantaneous_firing_rate"),
                    extensions=['.png', ])
        plt.close(f)
    except:
        traceback.print_exc()

    plot_analog_signal(all_exc_gsyn, variable_name="gsyn_exc",
                       ylabel="Exc synaptic conductance ($\mu S$)",
                       **common_values_for_plots)
    plot_analog_signal(all_inh_gsyn, variable_name="gsyn_inh",
                       ylabel="Inh synaptic conductance ($\mu S$)",
                       **common_values_for_plots)
    plot_analog_signal(all_voltages, variable_name="v",
                       ylabel="Membrane potential (mV)",
                       **common_values_for_plots)
    # THIS IS BROKEN!
    # print("Plotting a b c d post-synaptic hits")
    # l = 13
    # transparency_lvl = .7
    # for index, pop in enumerate(plot_order):
    #     if pop not in all_post_hits.keys():
    #         f = plt.figure(1, figsize=(l, l), dpi=400)
    #         plt.close(f)
    #         continue
    #     counts = all_post_hits[pop]
    #     f, axes = plt.subplots(2, 2, figsize=(l, l), dpi=400,
    #                            sharey=True, sharex=True)
    #     ax_a = axes[0, 0]
    #     ax_b = axes[0, 1]
    #     ax_c = axes[1, 0]
    #     ax_d = axes[1, 1]
    #
    #     ax_a.plot(counts.a, c="C0", alpha=transparency_lvl)
    #     ax_b.plot(counts.b, c="C1", alpha=transparency_lvl)
    #     ax_c.plot(counts.c, c="C2", alpha=transparency_lvl)
    #     ax_d.plot(counts.d, c="C3", alpha=transparency_lvl)
    #
    #     ax_a.set_title("a")
    #     ax_b.set_title("b")
    #     ax_c.set_title("c")
    #     ax_d.set_title("d")
    #
    #     ax_a.set_ylabel("# cases")
    #     ax_c.set_ylabel("# cases")
    #
    #     ax_c.set_ylabel("Time (ms)")
    #     ax_d.set_ylabel("Time (ms)")
    #
    #     plt.suptitle(use_display_name(pop))
    #
    #     plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
    #              stim_wanted_times.max() * time_to_bin_conversion)
    #     plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
    #     # plt.legend(loc="best")
    #     # plt.xlabel("Time (ms)")
    #
    #     plt.tight_layout()
    #     save_figure(
    #         plt,
    #         os.path.join(sim_fig_folder,
    #                      "post_hits_{}").format(pop),
    #         extensions=[".pdf", ])
    #     plt.close(f)

    print("Plotting a b c d post-synaptic hits V2")
    l = 13
    transparency_lvl = .7
    for index, pop in enumerate(plot_order):
        if pop not in all_post_hits.keys():
            f = plt.figure(1, figsize=(l, l), dpi=400)
            plt.close(f)
            continue

        sh_exc = np.max((other_recordings[pop]['v'].segments[0].filter(name='v')[0].magnitude.T * 2 ** 15).astype(int),
                        axis=0)
        sh_exc_2 = np.max(
            (other_recordings[pop]['gsyn_exc'].segments[0].filter(name='gsyn_exc')[0].magnitude.T * 2 ** 15).astype(
                int),
            axis=0)

        f = plt.figure(1, figsize=(l, l), dpi=400)
        a = ((sh_exc) & 0xFFFF)
        b = ((sh_exc >> 16) & 0xFFFF)
        c = ((sh_exc_2) & 0xFFFF)
        d = ((sh_exc_2 >> 16) & 0xFFFF)

        plt.plot(a + b + c + d, label='total', c="k", alpha=.3)
        plt.plot(d, label='d', c="C3", alpha=transparency_lvl)
        plt.plot(c, label='c', c="C2", alpha=transparency_lvl)
        plt.plot(b, label='b', c="C1", alpha=transparency_lvl)
        plt.plot(a, label='a', c="C0", alpha=transparency_lvl)

        plt.ylabel("Count")

        plt.suptitle(use_display_name(pop))

        plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
                 stim_wanted_times.max() * time_to_bin_conversion)
        plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
        plt.legend(loc="best")
        plt.xlabel("Time (ms)")

        plt.tight_layout()
        save_figure(
            plt,
            os.path.join(sim_fig_folder,
                         "or_post_hits_{}").format(pop),
            extensions=[".pdf", ])
        plt.close(f)

    print("Plotting worst_case spikes per pop (POTENTIALLY; ONLY IF TOOLS ARE SETUP UP THIS WAY")
    for _, pop in enumerate(plot_order):
        if (pop not in all_max_exc_spikes.keys() or
                pop not in all_max_inh_spikes.keys()):
            continue

        total = all_max_exc_spikes[pop] + all_max_inh_spikes[pop]
        f = plt.figure(1, figsize=(14, 9), dpi=400)
        plt.plot(total, c='k', label="Total", alpha=.3)

        plt.plot(all_max_exc_spikes[pop], rasterized=True,
                 label="excitatory spikes", alpha=transparency_lvl)
        plt.plot(all_max_inh_spikes[pop], rasterized=True,
                 label="inhibitory spikes", alpha=transparency_lvl)

        plt.title(use_display_name(pop))
        plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
                 stim_wanted_times.max() * time_to_bin_conversion)
        plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
        plt.legend(loc="best")
        plt.ylabel("Max # of spikes for a neuron")
        plt.xlabel("Time (ms)")
        plt.tight_layout()
        save_figure(
            plt,
            os.path.join(sim_fig_folder,
                         "potential_worst_neuron_only_spikes_per_pop_{}").format(pop),
            extensions=[".pdf", ])
        plt.close(f)

    if conn_exists and worst_case:
        for pop in plot_order:
            for proj, worst_spikes in per_conn_worst_spikes[pop].items():
                f = plt.figure(1, figsize=(12, 9), dpi=500)
                im = plt.imshow(worst_spikes,
                                interpolation='none',
                                extent=[stim_wanted_times.min() * time_to_bin_conversion,
                                        stim_wanted_times.max() * time_to_bin_conversion,
                                        0, worst_spikes.shape[1]],
                                origin='lower')
                ax = plt.gca()
                ax.set_aspect('auto')
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", "5%", pad="3%")
                cbar = plt.colorbar(im, cax=cax)
                cbar.set_label("# of spikes")

                ax.set_xlabel("Time (ms)")
                ax.set_ylabel("Neuron ID")

                ax.set_xlim(2800, 3800)
                #             plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
                #              stim_wanted_times.max() * time_to_bin_conversion)
                #             plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
                #             ax.xticks(np.linspace(2800, 3800, 6))
                ax.set_xticklabels([str(x) for x in (np.linspace(2800, 3800, 6) / 10)])

                save_figure(plt, os.path.join(sim_fig_folder,
                                              "afferent_wcs_{}_{}".format(pop, proj)),
                            #                         extensions=['.png', '.pdf'])
                            extensions=['.png', ])
                plt.show()
                plt.close(f)

    # Looking at number of spikes received per neuron for the purpose of computing Ring-Buffer left shifts to appropriately accomodate all (or most) spikes
    # Need to compute max, mean, 0.1 percentile and 99 percentile DURING STIMULUS PHASE -- Is this useful?
    spikes_per_proj_nutshell = {}
    if conn_exists and worst_case:
        excel_filename = os.path.join(sim_fig_folder,
                                      "per_projection_inc_spikes_description.xlsx")
        writer = pd.ExcelWriter(
            excel_filename,
            engine='xlsxwriter')
        for pop in plot_order:
            for proj, worst_spikes in per_conn_worst_spikes[pop].items():
                # need to filter worst spikes by timestep. Max wouldn't be affected, but all the other metrics would
                interval_of_interest = np.array(worst_spikes[:, 3000:3500])
                # remove zeros from this
                filterd_without_zeros = interval_of_interest[interval_of_interest != 0]
                percentile_01 = np.nanpercentile(filterd_without_zeros.ravel(), q=1)
                percentile_99 = np.nanpercentile(filterd_without_zeros.ravel(), q=99)
                mean = np.mean(filterd_without_zeros.ravel())
                std = np.std(filterd_without_zeros.ravel())
                maximums = np.max(filterd_without_zeros.ravel())
                spikes_per_proj_nutshell[proj] = [maximums, mean, std, percentile_01, percentile_99]

                # Also save a histogram
                f = plt.figure(1, figsize=(9, 9), dpi=400)
                plt.hist(filterd_without_zeros, bins=20, color=viridis_cmap(index / (n_plots + 1)),
                         rasterized=True,
                         edgecolor='k')

                plt.title(use_display_name(proj))

                plt.ylabel("Count")
                plt.xlabel("Spikes per timestep per neuron")
                plt.tight_layout()
                save_figure(
                    plt,
                    os.path.join(sim_fig_folder,
                                 "spike_count_hist_{}").format(proj),
                    extensions=[".pdf", ".png", ])
                plt.close(f)

        names = ['max', 'mean', 'std', '1st percentile', '99th percentile']
        df = pd.DataFrame(spikes_per_proj_nutshell, index=names)
        df.to_excel(writer, sheet_name="per_proj_inc_spikes_description")
        writer.save()

    # save required csv
    excel_filename = os.path.join(sim_fig_folder,
                                  "abcd_recordings_per_core.xlsx")

    writer = pd.ExcelWriter(
        excel_filename,
        engine='xlsxwriter')
    per_core_abcd = {}

    for index, pop in enumerate(plot_order):
        if pop not in per_conn_worst_spikes.keys():
            continue
        print("POPULATION", pop)
        zipper = zip((other_recordings[pop]['v'].segments[0].filter(name='v')[0].magnitude.T * 2 ** 15).astype(int),
                     (other_recordings[pop]['gsyn_exc'].segments[0].filter(name='gsyn_exc')[
                          0].magnitude.T * 2 ** 15).astype(int))

        abcd_data = []
        abcd_totals = []
        abcd_argmax = []
        for sh_exc, sh_exc_2 in zipper:
            x = np.asarray(unpack_abcd(sh_exc, sh_exc_2))
            total_x = np.sum(x, axis=0)
            max_x = np.argmax(total_x)
            abcd_data.append(x[:, max_x])
            abcd_totals.append(total_x)
            abcd_argmax.append(max_x)
        abcd_argmax = np.array(abcd_argmax)
        abcd_data = np.array(abcd_data)
        abcd_totals = np.array(abcd_totals)
        print("Number of cores", abcd_data.shape[0])
        names = ['a', 'b', 'c', 'd']  # N
        #     value = np.array(abcd_data, dtype=names)  # Note to self: numpy structured arrays are fucking broken
        array_for_pd = {name: abcd_data[:, col_id] for col_id, name in enumerate(names)}
        df = pd.DataFrame(array_for_pd)
        per_core_abcd[pop] = df
        df.to_excel(writer, sheet_name=pop)

    writer.save()

    transparency_lvl = .7
    for core_id, pop in enumerate(plot_order):
        if pop not in per_core_abcd.keys():
            #         f = plt.figure(1, figsize=(l, l), dpi=400)
            #         plt.close(f)
            continue
        f = plt.figure(1, figsize=(13, 10), dpi=400)
        abcd_df = per_core_abcd[pop]

        plt.plot(abcd_df[['a', 'b', 'c', 'd']].sum(axis=1), label='total',
                 c="k", ls="--", alpha=.6)
        plt.plot(abcd_df.d.values, label='d', c="C3", alpha=transparency_lvl)
        plt.plot(abcd_df.c.values, label='c', c="C2", alpha=transparency_lvl)
        plt.plot(abcd_df.b.values, label='b', c="C1", alpha=transparency_lvl)
        plt.plot(abcd_df.a.values, label='a', c="C0", alpha=transparency_lvl)

        plt.ylabel("Count")

        plt.title(use_display_name(pop))
        plt.legend(loc="best")
        plt.xlabel("Core")

        plt.tight_layout()
        save_figure(
            plt,
            os.path.join(sim_fig_folder,
                         "per_core_abcd_{}").format(pop),
            extensions=[".pdf", ])
        plt.show()
        plt.close(f)
    # Plot distribution of worst case spikes per population
    if conn_exists and worst_case:
        print("Plotting histogram of worst spike counts")
        for index, pop in enumerate(plot_order):
            counts = inc_spike_count[pop]
            maxs = np.max(counts, axis=1)
            assert maxs.size == all_neurons[pop]
            f = plt.figure(1, figsize=(9, 9), dpi=400)
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
                extensions=[".pdf", ])
            plt.close(f)

        print("Plotting worst_case spikes per pop")
        for _, pop in enumerate(plot_order):
            print("\t", pop)
            all_ws_contribution = per_conn_worst_spikes[pop]
            maxs = []
            labels = []
            for k, v in all_ws_contribution.items():
                labels.append(k)
                maxs.append(np.sum(v, axis=0))
            maxs = np.asarray(maxs)
            if len(maxs) == 0:
                continue
            total = np.sum(maxs, axis=0)
            assert total.size == maxs.shape[1]
            f = plt.figure(1, figsize=(14, 9), dpi=400)
            plt.plot(total, c='k', label="Total", alpha=.3)
            for i, row in enumerate(maxs):
                plt.plot(row, rasterized=True, label=labels[i], alpha=.7)
            plt.title(use_display_name(pop))

            plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
                     stim_wanted_times.max() * time_to_bin_conversion)
            plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
            plt.legend(loc="best")
            plt.ylabel("Max # of spikes per population")
            plt.xlabel("Time (ms)")
            plt.tight_layout()
            save_figure(
                plt,
                os.path.join(sim_fig_folder,
                             "sum_max_spikes_per_pop_{}").format(pop),
                extensions=[".pdf", ])
            plt.close(f)

            maxs = []
            labels = []
            for k, v in all_ws_contribution.items():
                labels.append(k)
                maxs.append(np.max(v, axis=0))
            maxs = np.asarray(maxs)
            if len(maxs) == 0:
                continue
            total = np.sum(maxs, axis=0)
            assert total.size == maxs.shape[1]
            f = plt.figure(1, figsize=(14, 9), dpi=400)
            plt.plot(total, c='k', label="Total", alpha=.3)
            for i, row in enumerate(maxs):
                plt.plot(row, rasterized=True, label=labels[i], alpha=.7)
            plt.title(use_display_name(pop))
            stim_wanted_times = np.linspace(time_filter[1] - 50,
                                            time_filter[2] + 50, 4).astype(int)
            plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
                     stim_wanted_times.max() * time_to_bin_conversion)
            plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
            plt.legend(loc="best")
            plt.ylabel("Max # of spikes for a neuron")
            plt.xlabel("Time (ms)")
            plt.tight_layout()
            save_figure(
                plt,
                os.path.join(sim_fig_folder,
                             "worst_neuron_only_spikes_per_pop_{}").format(pop),
                extensions=[".pdf", ])
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
    elephant_timestep_1 = sim_params_1['argparser']['timestep'] * pq.ms
    elephant_timestep_2 = sim_params_2['argparser']['timestep'] * pq.ms
    if elephant_timestep_1 != elephant_timestep_2:
        warnings.warn("The two simulation have different "
                      "simulation time steps dt",
                      UserWarning, stacklevel=2)
    ele_timesteps = [elephant_timestep_1, elephant_timestep_2]
    # TODO assert that parameters are the same between the simulations

    simtime = data_1['simtime'] * pq.ms
    useful_simtime = data_1['simtime'] * ms  # Use brian units as they are more usable
    assert (simtime == (data_2['simtime'] * pq.ms))

    same_sampling_rate = timestep == sim_params_2['argparser']['timestep'] * ms

    simulators = [sim_params_1['argparser']['simulator'],
                  sim_params_2['argparser']['simulator']]

    stimulus_params = data_1['stimulus_params'].ravel()[0]
    starts = np.cumsum(np.concatenate(([0], stimulus_params['stim_times']))) * pq.ms

    bin_size = 5 * pq.ms
    np_bin_size = 0.5 * pq.ms
    no_bins = int(simtime / bin_size)
    np_no_bins = int(simtime / np_bin_size)
    stimulus_no_bins = int((starts[2] - starts[1]) / bin_size)

    reporting_format_string = "{:30}-{:30}:{:6.3f}+-{:6.3f}"

    f = plt.figure(1, figsize=(9, 9), dpi=400)
    plt.close(f)

    all_spike_time_diff = {}
    all_spike_hist_diff = {}
    all_instantenous_rate_diff = {}
    all_corr_coef = {}
    all_cross_corr = {}
    all_covariances = {}
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
        cross_corr_per_neuron = []
        cov_per_neuron = []
        spike_time_tiling_coef_per_neuron = []
        van_rossum_dists = []
        spike_diff = []
        spike_hist_diff = []
        isis = {0: [],
                1: []}
        cvs = {0: [],
               1: []}
        cv2s = {0: [],
                1: []}

        # Comptue instantaneous rates for sim 1 and 2 and take the difference
        try:
            curr_inst_rates_1 = \
                elephant.statistics.instantaneous_rate(
                    pop_1_spikes,
                    sampling_period=elephant_timestep_1,
                    t_start=0 * pq.ms,
                    t_stop=simtime
                )

            curr_inst_rates_2 = \
                elephant.statistics.instantaneous_rate(
                    pop_2_spikes,
                    sampling_period=elephant_timestep_1,
                    t_start=0 * pq.ms,
                    t_stop=simtime
                )

            all_instantenous_rate_diff[pop] = np.expand_dims(np.squeeze(
                curr_inst_rates_1 - curr_inst_rates_2), 0)
        except ValueError as ve:
            traceback.print_exc()
            all_instantenous_rate_diff[pop] = [[np.nan]]

        # I guess we need to look at each neuron?
        for (p1, p2) in zip(pop_1_spikes, pop_2_spikes):
            # p1 = np.around(p1, 1) * pq.ms
            # p2 = np.around(p2, 1) * pq.ms
            try:
                spike_diff = p1 - p2
                # print("{:35} spike diff".format(pop), spike_diff)
            except:
                spike_diff = np.asarray([np.nan])

            np_bin_p1, np_bin_p1_edges = np.histogram(
                p1, bins=np.arange(np_no_bins) * np_bin_size)
            np_bin_p2, np_bin_p2_edges = np.histogram(
                p2, bins=np.arange(np_no_bins) * np_bin_size)

            spike_hist_diff.append(np_bin_p1 - np_bin_p2)

            # This section looks only at the stimulation period
            bin_p1 = elephant.conversion.BinnedSpikeTrain(
                p1, binsize=bin_size, num_bins=stimulus_no_bins,
                t_start=starts[1], t_stop=starts[2])
            bin_p2 = elephant.conversion.BinnedSpikeTrain(
                p2, binsize=bin_size, num_bins=stimulus_no_bins,
                t_start=starts[1], t_stop=starts[2])
            combined_binned = elephant.conversion.BinnedSpikeTrain(
                [p1, p2], binsize=bin_size, num_bins=stimulus_no_bins,
                t_start=starts[1], t_stop=starts[2])

            # This section looks at the entire simulation (from start to end
            full_bin_p1 = elephant.conversion.BinnedSpikeTrain(
                p1, binsize=bin_size, num_bins=no_bins,
                t_start=0 * pq.ms, t_stop=simtime)
            full_bin_p2 = elephant.conversion.BinnedSpikeTrain(
                p2, binsize=bin_size, num_bins=no_bins,
                t_start=0 * pq.ms, t_stop=simtime)
            full_combined_binned = elephant.conversion.BinnedSpikeTrain(
                [p1, p2], binsize=bin_size, num_bins=no_bins,
                t_start=0 * pq.ms, t_stop=simtime)

            cc_matrix = elephant.spike_train_correlation.corrcoef(
                combined_binned,
                fast=True)
            corr_coef_per_neuron.append(cc_matrix[0, 1])

            cov_matrix = elephant.spike_train_correlation.covariance(
                combined_binned,
                fast=True)
            cov_per_neuron.append(cov_matrix[0, 1])

            tiling = elephant.spike_train_correlation.spike_time_tiling_coefficient(
                p1, p2, dt=bin_size)
            spike_time_tiling_coef_per_neuron.append(tiling)

            dist = elephant.spike_train_dissimilarity.van_rossum_dist(
                [p1, p2], tau=bin_size)[0, 1]
            van_rossum_dists.append(dist)

            cross_corr = elephant.spike_train_correlation.cross_correlation_histogram(
                bin_p1, bin_p2,
                border_correction=False,
                binary=False, kernel=None, method='speed',
                cross_corr_coef=True)
            cross_corr_per_neuron.append(cross_corr)

            # Inter spike intervals
            isis[0].append(np.array(elephant.statistics.isi(
                p1 / pq.ms
            )))
            isis[1].append(np.array(elephant.statistics.isi(
                p2 / pq.ms
            )))

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
            "", "covariance",
            np.nanmean(cov_per_neuron),
            np.nanstd(cov_per_neuron)))

        print(reporting_format_string.format(
            "", "tiling coef",
            np.nanmean(spike_time_tiling_coef_per_neuron),
            np.nanstd(spike_time_tiling_coef_per_neuron)))

        print(reporting_format_string.format(
            "", "van rossum dist",
            np.nanmean(van_rossum_dists),
            np.nanstd(van_rossum_dists)))

        print(reporting_format_string.format(
            "", "instant. firing rate",
            np.nanmean(all_instantenous_rate_diff[pop]),
            np.nanstd(all_instantenous_rate_diff[pop])))

        all_corr_coef[pop] = corr_coef_per_neuron
        all_cross_corr[pop] = cross_corr_per_neuron
        all_covariances[pop] = cov_per_neuron
        all_distances[pop] = van_rossum_dists
        all_tiling_coef[pop] = spike_time_tiling_coef_per_neuron
        all_spike_time_diff[pop] = spike_diff
        all_spike_hist_diff[pop] = np.expand_dims(
            np.sum(spike_hist_diff, axis=0), 0)

        all_isi[pop] = isis
        all_cv[pop] = cvs
        all_cv2[pop] = cv2s

        # Log isis
        try:
            log_isi_1 = np.log(np.concatenate(all_isi[pop][0]))
            log_isi_2 = np.log(np.concatenate(all_isi[pop][1]))

            log_isi_1 = log_isi_1[np.isfinite(log_isi_1)]
            log_isi_2 = log_isi_2[np.isfinite(log_isi_2)]

            print("{:45}".format("log_isi_1 distribution normality"),
                  stats.normaltest(log_isi_1))
            print("{:45}".format("log_isi_2 distribution normality"),
                  stats.normaltest(log_isi_2))
            ttest_for_pop = stats.ttest_ind(log_isi_1, log_isi_2, equal_var=False)

            print("{:30}-{:30}:".format(
                "", "ttest on log ISI (pvalue, statistic)"),
                ttest_for_pop.pvalue, ttest_for_pop.statistic)
        except:
            traceback.print_exc()

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

    all_diff_error = {}
    if other_recordings_1 is not None and other_recordings_2 is not None:
        print("=" * 80)
        print("Analogue value analysis (avg +- std)")
        # https://elephant.readthedocs.io/en/latest/reference/spectral.html
        print("-" * 80)
        for signal in ['v', 'gsyn_exc', 'gsyn_inh']:
            all_coherence[signal] = {}
            all_freqs[signal] = {}
            all_lags[signal] = {}

            all_diff_error[signal] = {}
            for pop in plot_order:
                if pop not in other_recordings_1.keys():
                    continue
                pop_1_v = other_recordings_1[pop][signal].segments[0].analogsignals
                pop_2_v = other_recordings_2[pop][signal].segments[0].analogsignals

                welch_cohere_freqs_per_neuron = []
                welch_cohere_lags_per_neuron = []
                welch_cohere_coher_per_neuron = []

                diff_per_neuron = []
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

                    if same_sampling_rate:
                        diff_per_neuron.append(p1 - p2)

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
                all_lags[signal][pop] = welch_cohere_lags_per_neuron
                all_freqs[signal][pop] = welch_cohere_freqs_per_neuron

                if same_sampling_rate:
                    all_diff_error[signal][pop] = np.asarray(diff_per_neuron)

    print("=" * 80)
    print("Plotting figures...")
    print("-" * 80)

    # raster plot including ALL populations
    print("Plotting spiking raster plot for all populations")

    # Pre-compute conversions
    time_to_bin_conversion = 1. / (timestep / ms)
    bins_in_3ms = int(3 * time_to_bin_conversion)
    no_timesteps = int(simtime / ms * time_to_bin_conversion)
    wanted_times = np.linspace(0, (useful_simtime / ms), 6).astype(int)

    common_highlight_values = {
        'start': None,  # stim_period_start,
        'stop': None,  # stim_period_end,
        'increment': timestep / ms,
    }

    common_values_for_plots = {
        'plot_order': plot_order,
        'wanted_times': wanted_times,  # wanted_times,
        'time_to_bin_conversion': time_to_bin_conversion,  # time_to_bin_conversion,
        'fig_folder': sim_fig_folder,
        'highlight_stim': None,  # highlight_stim,
        'common_highlight_values': None,  # common_highlight_values,
    }

    plot_analog_signal(all_instantenous_rate_diff,
                       "instantaneous_rate_diff",
                       ylabel="Difference (Hz)",
                       **common_values_for_plots)

    plot_analog_signal(all_spike_hist_diff,
                       "spike_hist_diff",
                       ylabel="Difference (Count)",
                       scale_xticks=1 / (10 * np_bin_size / pq.ms),
                       **common_values_for_plots)

    if same_sampling_rate:
        plot_analog_signal(all_diff_error['v'],
                           "v_diff_error",
                           ylabel="Error (mV)",
                           **common_values_for_plots)
        plot_analog_signal(all_diff_error['gsyn_exc'],
                           "gsyn_exc_diff_error",
                           ylabel="Error ($\mu S$)",
                           **common_values_for_plots)
        plot_analog_signal(all_diff_error['gsyn_inh'],
                           "gsyn_inh_diff_error",
                           ylabel="Error ($\mu S$)",
                           **common_values_for_plots)

    # Plot various analog variables
    for signal in all_coherence.keys():
        plot_analog_signal_with_error(all_coherence[signal], signal + "_coherence",
                                      xlabel="Frequency (Hz)",
                                      ylabel="Coherence",
                                      xticks=all_freqs[signal],
                                      errorbars=False,
                                      **common_values_for_plots)

        plot_analog_signal_with_error(all_lags[signal], signal + "_lag",
                                      xlabel="Frequency (Hz)",
                                      ylabel="Lag (ms)",
                                      xticks=all_freqs[signal],
                                      errorbars=False,
                                      **common_values_for_plots)
    # Plotting boxplots
    print("PLOTTING", "all_corr_coef")
    plot_boxplot_for_all_pops(all_corr_coef, "corr_coeff",
                              xlabel="Population",
                              ylabel="Correlation Coefficient",
                              **common_values_for_plots)

    try:

        plot_histogram_for_all_pops(all_cross_corr, "cross_corr",
                                    xlabel="Population",
                                    ylabel="Cross Correlation",
                                    **common_values_for_plots)
        plot_boxplot_for_all_pops(all_cross_corr, "cross_corr",
                                  xlabel="Population",
                                  ylabel="Cross Correlation",
                                  **common_values_for_plots)
    except:
        pass

    plot_boxplot_for_all_pops(all_covariances, "covariance",
                              xlabel="Population",
                              ylabel="Covariance",
                              **common_values_for_plots)

    plot_boxplot_for_all_pops(all_distances, "van_rossum_dist",
                              xlabel="Population",
                              ylabel="van Rossum Distance",
                              **common_values_for_plots)

    plot_boxplot_for_all_pops(all_tiling_coef, "tiling_coeff",
                              xlabel="Population",
                              ylabel="Tiling Coefficient",
                              **common_values_for_plots)
    # Plot side by side boxplots
    plot_sbs_boxplot_for_all_pops(all_isi, "isi",
                                  xlabel="Population",
                                  ylabel="ISI (ms)", titles=simulators,
                                  yscale='log',
                                  **common_values_for_plots)

    plot_sbs_boxplot_for_all_pops(all_cv, "cv",
                                  xlabel="Population",
                                  ylabel="CV", titles=simulators,
                                  **common_values_for_plots)

    plot_sbs_boxplot_for_all_pops(all_cv2, "cv2",
                                  xlabel="Population",
                                  ylabel="CV2", titles=simulators,
                                  **common_values_for_plots)

    # Plot side by side histograms
    # plot_sbs_histogram_for_all_pops(all_isi, "isi",
    #                                 xlabel="ISI",
    #                                 ylabel="count", titles=simulators,
    #                                 **common_values_for_plots)
    #
    # plot_sbs_histogram_for_all_pops(all_cv, "cv",
    #                                 xlabel="CV",
    #                                 ylabel="count", titles=simulators,
    #                                 **common_values_for_plots)
    #
    # plot_sbs_histogram_for_all_pops(all_cv2, "cv2",
    #                                 xlabel="CV2",
    #                                 ylabel="count", titles=simulators,
    #                                 **common_values_for_plots)

    # Plotting histograms
    # plot_histogram_for_all_pops(all_cross_corr, "cross_corr",
    #                             xlabel="Cross Correlation",
    #                             ylabel="count",
    #                             **common_values_for_plots)

    # plot_histogram_for_all_pops(all_corr_coef, "corr_coeff",
    #                             xlabel="Correlation Coefficient",
    #                             ylabel="count",
    #                             **common_values_for_plots)
    #
    # plot_histogram_for_all_pops(all_distances, "van_rostum_dist",
    #                             xlabel="van Rostum Distance",
    #                             ylabel="count",
    #                             **common_values_for_plots)
    #
    # plot_histogram_for_all_pops(all_tiling_coef, "tiling_coeff",
    #                             xlabel="Tiling Coefficient",
    #                             ylabel="count",
    #                             **common_values_for_plots)

    print("=" * 80)


def scatter_plot(data, variable_name,
                 xlabel, ylabel, plot_order,
                 wanted_times, time_to_bin_conversion,
                 fig_folder,
                 highlight_stim, common_highlight_values,
                 titles):
    print("Scatter plot -- {} for all populations".format(
        variable_name))
    f, axes = plt.subplots(len(data.keys()), 1,
                           figsize=(14, 20), sharex=True, dpi=400)
    n_plots = len(plot_order)
    minimus = 0
    maximus = 0
    for index, pop in enumerate(plot_order):
        curr_ax = axes[index]
        curr_pop = data[pop]
        curr_ax.scatter(np.asarray(curr_pop[0]).ravel(),
                        np.asarray(curr_pop[1]).ravel(),
                        color=color_for_index(index, n_plots),
                        rasterized=True)
        if index == 0:
            curr_ax.set_title(use_display_name(titles[0]))
        curr_ax.set_ylabel(use_display_name(pop) + "\n" + ylabel)
        curr_ax.grid(True, which="major", axis="both")

        minimus = min(minimus, min(data[pop][0]), min(data[pop][1]))
        maximus = max(maximus, max(data[pop][0]), max(data[pop][1]))

    plt.xlabel(xlabel)
    plt.xlim([minimus - (.1 * maximus), 1.1 * maximus])
    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "scatter_plot_{}".format(variable_name)),
                extensions=['.png', '.pdf'])
    plt.close(f)


def plot_analog_signal_with_error(data, variable_name,
                                  xlabel, ylabel, plot_order,
                                  wanted_times, time_to_bin_conversion,
                                  fig_folder,
                                  highlight_stim, common_highlight_values,
                                  errorbars=True,
                                  xticks=None):
    print("Plotting {} traces for each population with errors".format(variable_name))
    n_plots = len(plot_order)
    for index, pop in enumerate(plot_order):
        if pop in ["glomerulus"] or pop not in data.keys():
            f = plt.figure(1, figsize=(9, 9), dpi=400)
            plt.close(f)
            continue

        values_for_pop = np.squeeze(data[pop])
        if len(values_for_pop.shape) == 1:
            values_for_pop = np.expand_dims(values_for_pop, 1)
        if xticks is not None:
            adjusted_xticks = xticks[pop].ravel()
            argsorted_ticks = np.argsort(
                adjusted_xticks[np.logical_and(adjusted_xticks > 0, adjusted_xticks < 2000)])
            adjusted_xticks = adjusted_xticks[argsorted_ticks]
            passed_in_ticks = True
        else:
            adjusted_xticks = np.arange(values_for_pop.shape[1])
            passed_in_ticks = False
            argsorted_ticks = adjusted_xticks
        f = plt.figure(1, figsize=(9, 9), dpi=400)
        if errorbars:
            plt.errorbar(adjusted_xticks,
                         np.nanmean(values_for_pop, axis=1).ravel()[argsorted_ticks],
                         yerr=np.nanstd(values_for_pop, axis=1).ravel(),
                         color=color_for_index(index, n_plots),
                         rasterized=True)
        else:
            plt.plot(adjusted_xticks,
                     np.nanmean(values_for_pop, axis=1).ravel()[argsorted_ticks],
                     color=color_for_index(index, n_plots),
                     rasterized=True)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        save_figure(plt, os.path.join(fig_folder,
                                      "{}_{}".format(pop, variable_name)),
                    extensions=['.png', '.pdf'])
        plt.close(f)

        # Also plot using matshow
        if not passed_in_ticks:
            extent_xticks = [np.min(wanted_times), np.max(wanted_times)]
        else:
            extent_xticks = [np.min(adjusted_xticks), np.max(adjusted_xticks)]

        # values_for_pop = np.squeeze(values_for_pop)
        f = plt.figure(1, figsize=(18, 9), dpi=500)
        im = plt.imshow(values_for_pop.T,
                        interpolation='none',
                        extent=[extent_xticks[0], extent_xticks[1],
                                0, values_for_pop.shape[1]],
                        origin='lower')
        ax = plt.gca()
        ax.set_aspect('auto')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(ylabel)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Neuron ID")
        save_figure(plt, os.path.join(fig_folder,
                                      "imshow_{}_{}".format(pop, variable_name)),
                    extensions=['.png', '.pdf'])
        plt.close(f)


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
    n_plots = len(plot_order)

    bp_width = 0.7
    f = plt.figure(figsize=(12, 8), dpi=600)
    for index, pop in enumerate(plot_order):
        curr_data = np.asarray(data[pop])
        curr_data = curr_data[np.isfinite(curr_data)]
        plt.boxplot(curr_data, notch=True, positions=[index + 1],
                    medianprops=dict(
                        color=color_for_index(index, n_plots),
                        linewidth=1.5),
                    widths=bp_width
                    )
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
                                  titles, yscale='linear'):
    print("Plotting SIDE BY SIDE boxplot for {} for all populations".format(
        variable_name))
    n_plots = len(plot_order)
    spacer = 4
    bp_width = 0.7
    f = plt.figure(figsize=(12, 8), dpi=600)
    for index, pop in enumerate(plot_order):
        curr_pop = data[pop]
        curr_res_0 = np.asarray(curr_pop[0])
        curr_res_1 = np.asarray(curr_pop[1])
        # Filter NaNs
        curr_res_0 = curr_res_0[np.isfinite(curr_res_0)]
        curr_res_1 = curr_res_1[np.isfinite(curr_res_1)]
        bp0 = plt.boxplot(curr_res_0.ravel(), notch=True,
                          positions=[spacer * index + 1],
                          medianprops=dict(color=viridis_cmap(0.0),
                                           linewidth=1.5),
                          widths=bp_width)
        bp1 = plt.boxplot(curr_res_1.ravel(), notch=True,
                          positions=[spacer * index + 2],
                          medianprops=dict(color=viridis_cmap(.7),
                                           linewidth=1.5),
                          widths=bp_width)
    title_display_names = [use_display_name(x) for x in titles]
    plt.legend([bp0["medians"][0], bp1["medians"][0]], title_display_names,
               loc='best')
    plt.ylabel(ylabel)
    plt.xlim([0, spacer * n_plots])
    plt.grid(True, which="major", axis="y")
    plt.xlabel(xlabel)
    plt.yscale(yscale)
    if yscale != 'linear':
        suffix = "_" + yscale
    else:
        suffix = ""
    xtick_display_names = [use_display_name(x) for x in plot_order]
    _, labels = plt.xticks(spacer * np.arange(len(plot_order)) + 1.5,
                           xtick_display_names)

    f.tight_layout()
    save_figure(plt, os.path.join(fig_folder,
                                  "bp_sbs_all_pop_{}".format(variable_name + suffix)),
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
