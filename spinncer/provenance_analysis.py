from spinncer.analysis_common import *
from os.path import join as join


def extract_per_pop_placements(df, pops):
    placement_results = {}
    for pop in pops:
        pop_df = df[df['pop'] == pop]
        placement_df = pop_df[["x", "y", "p", "no_atoms"]].drop_duplicates()
        placement_results[pop] = placement_df
    return placement_results


def extract_per_pop_info(df, type_of_prov, pops, report=False):
    pop_results = {k: None for k in pops}
    pop_results['global_mean'] = np.nan
    pop_results['global_max'] = np.nan
    pop_results['global_min'] = np.nan

    prov_filter = df['prov_name'] == type_of_prov
    filtered_prov = df[prov_filter]
    if report:
        print("{:40} for populations:".format(type_of_prov))

    _means = []
    _maxs = []
    _mins = []
    for pop in pops:
        pop_df = filtered_prov[filtered_prov['pop'] == pop]
        curr_pop_values = pop_df.prov_value
        _mean = curr_pop_values.mean()
        _max = curr_pop_values.max()
        _min = curr_pop_values.min()
        if report:
            print("\t{:25} - avg {:10.2f} max {:10.2f}".format(
                pop, curr_pop_values.mean(), curr_pop_values.max()))
        # save values
        pop_results[pop] = {
            'mean': _mean,
            'max': _max,
            'min': _min,
            'all': curr_pop_values
        }
        _means.append(_mean)
        _maxs.append(_max)
        _mins.append(_min)
    if report:
        write_line()
    pop_results['global_mean'] = np.nanmean(np.asarray(_means))
    pop_results['global_max'] = np.nanmax(np.asarray(_maxs))
    pop_results['global_min'] = np.nanmin(np.asarray(_mins))
    return pop_results


def provenance_csv_analysis(in_folder, fig_folder):
    write_header("Reading provenances in folder " + in_folder)
    prov = pd.read_csv(join(in_folder, "structured_provenance.csv"))
    pops = prov['pop'].unique()
    pops.sort()
    types_of_provenance = prov['prov_name'].unique()
    prov_of_interest = [
        'MAX_SPIKES_IN_A_TICK',
        'Times_synaptic_weights_have_saturated',
        'late_packets',
        'Times_the_input_buffer_lost_packets',
        'Times_the_timer_tic_over_ran'
    ]

    results = {k: None for k in types_of_provenance}
    placements = {}
    # TODO report number of neurons to make sure the networks is correct
    write_short_msg("DETECTED POPULATIONS", pops)

    for type_of_prov in types_of_provenance:
        rep = True if type_of_prov in prov_of_interest else False
        results[type_of_prov] = extract_per_pop_info(prov, type_of_prov, pops,
                                                     report=rep)
    placements = extract_per_pop_placements(prov, pops)
    return results, types_of_provenance, prov_of_interest, placements


def sweep_provenance_analysis(in_folder, fig_folder, group_on):
    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)
    current_fig_folder = join(fig_folder, in_folder)
    # Make folder for current figures
    if not os.path.isdir(current_fig_folder) and not os.path.exists(current_fig_folder):
        os.mkdir(current_fig_folder)

    write_header("Analysing provenances in folder " + in_folder)
    folder_contents = os.listdir(in_folder)
    batch_infos = [fc for fc in folder_contents if "batch_" in fc]
    if len(batch_infos) > 1:
        raise FileExistsError("There are multiple batch info .npz archives "
                              "in the current directory!")
    if len(batch_infos) == 1:
        batch_info = np.load(join(in_folder, batch_infos[0]))
        run_id = "@" + batch_infos[0].split("_")[1].split(".")[0]
        if batch_info:
            write_header("Batch info contained in " + batch_infos[0])
            write_short_msg("FILES", batch_info.files)
            poi = batch_info['parameters_of_interest'].ravel()[0]
            calls = batch_info['log_calls']
            write_short_msg("PARAMETERS", poi)
        run_folders = [fc for fc in folder_contents if run_id in fc]

    if len(batch_infos) == 0:
        print("No batch info detected. Going to assume that "
              "all folders in the current directory are "
              "related to the run.")
        calls, poi = None
        run_folders = [fc for fc in folder_contents if os.path.isdir(join(in_folder, fc))]
    collated_results = {k: None for k in run_folders}
    types_of_provenance = None
    prov_of_interest = None
    placements = {}
    write_short_msg("Number of folders", len(run_folders))
    write_sep()
    for folder in run_folders:
        collated_results[folder], types_of_provenance, \
        prov_of_interest, placements[folder] = provenance_csv_analysis(
            join(in_folder, folder), fig_folder)

    if group_on is None:
        write_header("REPORTING BEST SIMULATIONS")
        cumulative_report(collated_results, types_of_provenance,
                          prov_of_interest)
    else:
        for condition in group_on:
            write_header("GROUPING ON CONDITION {}".format(condition))
            new_collated_results = {key: value for (key, value) in
                                    collated_results.items() if condition in key}
            cumulative_report(new_collated_results, types_of_provenance,
                              prov_of_interest)

    plot_population_placement(collated_results, placements,
                              fig_folder=current_fig_folder)

    plot_per_population_max_spikes_per_tick(collated_results, calls, poi,
                                            prov_of_interest,
                                            group_on, current_fig_folder)


def cumulative_report(collated_results, types_of_provenance, prov_of_interest):
    sorted_key_list = list(collated_results.keys())
    sorted_key_list.sort()
    for type_of_prov in types_of_provenance:
        of_interest = True if type_of_prov in prov_of_interest else False
        if not of_interest:
            continue
        # Need to loop over all dicts in collated_results
        _means = []
        _maxs = []
        _mins = []
        _keys = []
        _max_values_per_pop = None

        print("{:40} for all cases".format(type_of_prov))
        for k in sorted_key_list:
            v = collated_results[k]
            filtered_v = v[type_of_prov]
            _keys.append(k)
            _means.append(filtered_v['global_mean'])
            _maxs.append(filtered_v['global_max'])
            _mins.append(filtered_v['global_min'])
            print("{:40} | min {:10.2f} | mean {:10.2f} | max {:10.2f}".format(
                k.split("@")[0], filtered_v['global_min'], filtered_v['global_mean'], filtered_v['global_max']
            ))
            if _max_values_per_pop is None:
                _max_values_per_pop = {x: [] for x in filtered_v.keys() if "cell" in x}
            for vp in _max_values_per_pop.keys():
                _max_values_per_pop[vp].append(filtered_v[vp]['max'])
        # Report cumulative stats per population
        write_line()
        print("{:40} for all populations".format(type_of_prov))
        reporting_keys = list(_max_values_per_pop.keys())
        reporting_keys.sort()
        for rk in reporting_keys:
            vals = _max_values_per_pop[rk]
            print("{:40} | mean {:10.2f} | max {:10.2f} | std {:10.2f}".format(
                rk, np.nanmean(vals), np.nanmax(vals), np.nanstd(vals)
            ))


def plot_population_placement(collated_results, placements, fig_folder):
    write_header("PLOTTING MAPS")
    sorted_key_list = list(collated_results.keys())
    sorted_key_list.sort()
    for selected_sim in sorted_key_list:
        filtered_placement = \
            placements[selected_sim]
        placements_per_pop = {x: filtered_placement[x]
                              for x in filtered_placement.keys()
                              if "cell" in x}

        # Plotting bit
        # Fake printing to start things off...
        f = plt.figure(1, figsize=(9, 9), dpi=400)
        plt.close(f)
        # Compute plot order
        plot_order = get_plot_order(placements_per_pop.keys())
        plot_display_names = []
        for po in plot_order:
            plot_display_names.append(use_display_name(po))
        n_plots = len(plot_order)

        collated_placements = pd.concat([
            filtered_placement[x] for x in plot_order
        ])

        max_x = collated_placements.x.max() * 5
        max_y = collated_placements.y.max() * 5
        x_ticks = np.arange(0, max_x, 5)[::2]
        # x_tick_lables = np.linspace(0, collated_placements.x.max(), 6).astype(int)
        x_tick_lables = (x_ticks / 5).astype(int)
        y_ticks = np.arange(0, max_y, 5)[::2]
        # y_tick_lables = np.linspace(0, collated_placements.y.max(), 6).astype(int)
        y_tick_lables = (y_ticks / 5).astype(int)
        map = np.ones((max_x, max_y)) * np.nan
        for index, pop in enumerate(plot_order):
            curr_pl = placements_per_pop[pop]
            for row_index, row in curr_pl.iterrows():
                map[
                    int(4 * row.x + (row.p // 4)),
                    int(4 * row.y + (row.p % 4))
                ] = index

        uniques = np.unique(map[np.isfinite(map)]).astype(int)
        # norm = mlib.colors.BoundaryNorm(uniques, n_plots)
        # fmt = mlib.ticker.FuncFormatter(lambda x, pos: plot_order[::-1][norm(x)])
        # crop_point = np.max(np.max(np.argwhere(np.isfinite(map)), axis=0))
        f = plt.figure(1, figsize=(9, 9), dpi=500)
        # plt.matshow(map[:crop_point, :crop_point], interpolation='none')
        im = plt.imshow(map, interpolation='none', vmin=0, vmax=n_plots,
                        cmap=plt.get_cmap('viridis', n_plots))
        ax = plt.gca()

        plt.xticks(x_ticks, x_tick_lables)
        plt.yticks(y_ticks, y_tick_lables)
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(5))


        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label("Population")
        cbar.ax.set_yticks(uniques)
        cbar.ax.set_yticklabels(plot_display_names)

        save_figure(plt, join(fig_folder,
                              "map_of_placements_for_{}".format(selected_sim)),
                    extensions=['.png', '.pdf'])
        plt.close(f)

        # Some reports
        write_short_msg("Plotting map for", selected_sim)
        write_short_msg("Number of cores used", collated_placements.shape[0])
        write_short_msg("Number of chips used",
                        collated_placements[["x", "y"]].drop_duplicates().shape[0])
        write_short_msg("Unique pop ids", uniques)
        write_line()

def plot_per_population_max_spikes_per_tick(collated_results, calls, poi, prov,
                                            group_on, fig_folder):
    write_header("PLOTTING PER POPULATION VALUES FOR PROVENANCE")
    sorted_key_list = list(collated_results.keys())
    sorted_key_list.sort()
    f = plt.figure(1, figsize=(9, 9), dpi=400)
    plt.close(f)
    # Look at parameters in calls
    # Check if group_on is a parameter in poi and calls
    # Grab all of the uniques from
    for curr_g in group_on:
        # if curr_g is in POI then we can proceed with plotting things like this
        if curr_g in poi.keys():
            curr_poi = poi[curr_g]
            print("Current parameter of interest is ", curr_g,
                  "with values", curr_poi)
            # sort, just in case
            curr_poi.sort()
            for type_of_prov in prov:
                # Now we need to start collecting values
                # match curr_poi value with collated_results.keys() via calls
                curr_mapping = {}
                for curr_poi_val in curr_poi:
                    for call in calls:
                        if call[2][curr_g] == curr_poi_val:
                            # store filename associated with current value
                            match_fname = call[1]
                            filtered_collated_results = \
                                collated_results[match_fname][type_of_prov]
                            _max_values_per_pop = None
                            if _max_values_per_pop is None:
                                _max_values_per_pop = {x: None for x in filtered_collated_results.keys() if "cell" in x}
                            for vp in _max_values_per_pop.keys():
                                _max_values_per_pop[vp] = filtered_collated_results[vp]['max']
                            curr_mapping[curr_poi_val] = _max_values_per_pop
                            break
                # Plotting bit
                plot_order = get_plot_order(_max_values_per_pop.keys())
                n_plots = float(len(plot_order))

                f = plt.figure(1, figsize=(9, 9), dpi=400)
                for index, pop in enumerate(plot_order):
                    vals = []
                    for k in curr_poi:
                        vals.append(curr_mapping[k][pop])
                    # vals = np.sort(np.asarray(vals).view('i8,i8'), order=['f1'], axis=0).view(np.int)
                    vals = np.array(vals)
                    if np.any(np.isfinite(vals)):
                        plt.plot(curr_poi, vals,
                                 color=color_for_index(index, n_plots),
                                 marker='o',
                                 label=use_display_name(pop),
                                 alpha=0.8)

                plt.xlabel(use_display_name(curr_g))
                plt.ylabel(use_display_name(type_of_prov))
                plt.legend(loc='best')
                plt.tight_layout()
                save_figure(plt, join(fig_folder, "{}".format(type_of_prov)),
                            extensions=['.png', '.pdf'])
                plt.close(f)


if __name__ == "__main__":
    from spinncer.analysis_argparser import *

    if analysis_args.input and len(analysis_args.input) > 0:
        for in_folder in analysis_args.input:
            sweep_provenance_analysis(in_folder, analysis_args.figures_dir,
                                      group_on=analysis_args.group_on)
    else:
        # Constants
        fig_folder = "figures"
        pass
