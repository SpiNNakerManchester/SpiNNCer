import itertools

import pandas as pd
from spinncer.analysis_common import *
from os.path import join as join
from numpy.polynomial.polynomial import Polynomial
import traceback

from spinn_front_end_common.utilities.globals_variables import (
    provenance_file_path)
from spinn_front_end_common.interface.provenance import ProvenanceReader


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


def save_provenance_to_file_from_database(in_file, simulator):
    # Here we need to get the provenance from the database and put it in
    # the specified file

    # list provenance of interest
    router_provenance_of_interest = [
        'Dumped_from_a_Link',
        'Dumped_from_a_processor',
        'Local_Multicast_Packets',
        'External_Multicast_Packets',
        'Dropped_Multicast_Packets',
        'Missed_For_Reinjection'
    ]
    prov_of_interest = [
        'Maximum number of spikes in a timer tick',
        'Times_synaptic_weights_have_saturated',
        'late_packets',
        'Times_the_input_buffer_lost_packets',
        'Times_the_timer_tic_over_ran',
        'Total_pre_synaptic_events',
        'Maximum number of DMAs in a timer tick',
        'Maximum pipeline restarts',
        'send_multicast_packets',
        'Maximum number of spikes flushed in a timer tick',
        'Total number of spikes flushed'
    ]

    # Custom provenance presentation from SpiNNCer
    # write provenance to file here in a useful way
    columns = ['pop', 'label', 'min_atom', 'max_atom', 'no_atoms',
               'x', 'y', 'p',
               'prov_name', 'prov_value',
               'fixed_sdram', 'sdram_per_timestep',
               'cpu_cycles', 'dtcm']
    # assert (len(prov_placement) == len(prov_items))
    structured_provenance = list()
    metadata = {}
    # Retrieve filename from spynnaker8/spinnaker.py
    provenance_filename = in_file

    if provenance_filename:
        # Produce metadata from the simulator info
        metadata['name'] = simulator.name
        metadata['no_machine_time_steps'] = simulator.no_machine_time_steps
        metadata['machine_time_step'] = simulator.machine_time_step
        # metadata['config'] = simulator.config
        metadata['machine'] = simulator.machine
        metadata['structured_provenance_filename'] = \
            simulator.structured_provenance_filename

        # how do we loop over all the placements in the database at this point?
        # can we get router_provenance from the database?
        # placements is 0... so is this router(0, 0)?
        pr = ProvenanceReader(os.path.join(
            provenance_file_path(), "provenance.sqlite3"))

        cores_list = pr.get_cores_with_provenance()

        # router_provenance = pr.get_provenance_for_router(0, 0)
        # print("router_provenance: ", router_provenance)

        # for i, (provenance, placement) in enumerate(zip(prov_items, prov_placement)):
        for core in cores_list:
            x = core[1]
            y = core[2]
            p = core[3]
            structured_prov_core = get_provenance_for_core(pr, x, y, p)

            pop = structured_prov_core['pop']
            if pop == []:
                continue
            pop = pop[0][0]
            fixed_sdram = structured_prov_core['fixed_sdram'][0][0]
            sdram_per_timestep = structured_prov_core['sdram_per_timestep'][0][0]
            cpu_cycles = structured_prov_core['cpu_cycles'][0][0]

            label = structured_prov_core['label'][0][0]
            max_atom = structured_prov_core['max_atom'][0][0]
            min_atom = structured_prov_core['min_atom'][0][0]
            no_atoms = structured_prov_core['no_atoms'][0][0]
            dtcm = structured_prov_core['dtcm'][0][0]

            for prov_name in prov_of_interest:
                prov_value = get_core_provenance_value(pr, x, y, p, prov_name)
                if prov_value == []:
                    prov_value = 0
                else:
                    prov_value = prov_value[0][0]

                structured_provenance.append(
                    [pop, label, min_atom, max_atom, no_atoms,
                     x, y, p,
                     prov_name, prov_value,
                     fixed_sdram, sdram_per_timestep,
                     cpu_cycles, dtcm]
                )

            for prov_name in router_provenance_of_interest:
                prov_value = get_router_provenance_value(pr, x, y, prov_name)
                if prov_value == []:
                    prov_value = 0
                else:
                    prov_value = prov_value[0][0]

                structured_provenance.append(
                    [pop, label, min_atom, max_atom, no_atoms,
                     x, y, p,
                     prov_name, prov_value,
                     fixed_sdram, sdram_per_timestep,
                     cpu_cycles, dtcm]
                )

        # print("structured provenance: ", structured_provenance)

        structured_provenance_df = pd.DataFrame.from_records(
            structured_provenance, columns=columns)

        structured_provenance_df.to_csv("structured_provenance.csv")


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
        'Times_the_timer_tic_over_ran',
        'Total_pre_synaptic_events',
        'MAX_DMAS_IN_A_TICK',
        'MAX_PIPELINE_RESTARTS',
        'send_multicast_packets',
        'MAX_FLUSHED_SPIKES',
        'TOTAL_FLUSHED_SPIKES'
    ]

    results = {k: None for k in types_of_provenance}
    # TODO report number of neurons to make sure the networks is correct
    write_short_msg("DETECTED POPULATIONS", pops)

    for type_of_prov in types_of_provenance:
        rep = True if type_of_prov in prov_of_interest else False
        results[type_of_prov] = extract_per_pop_info(prov, type_of_prov, pops,
                                                     report=rep)
    placements = extract_per_pop_placements(prov, pops)
    return results, types_of_provenance, prov_of_interest, placements


def get_provenance_for_core(pr, x, y, p):
    structured_prov = {}
    columns_to_get = ['pop', 'label', 'min_atom', 'max_atom', 'no_atoms',
               'fixed_sdram', 'sdram_per_timestep',
               'cpu_cycles', 'dtcm']  # add more as needed

    for column_to_get in columns_to_get:
        query = """
            SELECT the_value
            FROM core_provenance
            WHERE x = ? AND y = ? AND p = ? AND description = ?
            """
        # result = pr.run_query(query, [x, y, p, column_to_get])
        structured_prov[column_to_get] = pr.run_query(
            query, [x, y, p, column_to_get])
        # print('x,y,p,desc: ', x, y, p, column_to_get,
        #       structured_prov[column_to_get])

    return structured_prov


def get_core_provenance_value(pr, x, y, p, description):
    query = """
        SELECT the_value
        FROM core_provenance
        WHERE x = ? AND y = ? AND p = ? AND description = ?
        """
    return pr.run_query(query, [x, y, p, description])


def get_router_provenance_value(pr, x, y, description):
    query = """
        SELECT the_value
        FROM router_provenance
        WHERE x = ? AND y = ? AND description = ?
        """
    return pr.run_query(query, [x, y, description])


def sweep_provenance_analysis(in_folder, fig_folder, group_on,
                              group_on_name):
    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)
    current_fig_folder = join(fig_folder, in_folder.split('/')[-1])
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
        batch_info = np.load(join(in_folder, batch_infos[0]), allow_pickle=True)
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

    router_provenance_of_interest = [
        'Dumped_from_a_Link',
        'Dumped_from_a_processor',
        'Local_Multicast_Packets',
        'External_Multicast_Packets',
        'Dropped_Multicast_Packets',
        'Missed_For_Reinjection'
    ]
    router_pop_names = ['router_provenance']

    write_short_msg("Number of folders", len(run_folders))
    write_sep()
    for folder in run_folders:
        collated_results[folder], types_of_provenance, \
        prov_of_interest, placements[folder] = provenance_csv_analysis(
            join(in_folder, folder), fig_folder)

    if group_on_name is None:
        write_header("REPORTING BEST SIMULATIONS")
        cumulative_report(collated_results, types_of_provenance,
                          prov_of_interest)
    else:
        for condition in group_on_name:
            write_header("GROUPING ON CONDITION {}".format(condition))
            new_collated_results = {key: value for (key, value) in
                                    collated_results.items() if condition in key}
            cumulative_report(new_collated_results, types_of_provenance,
                              prov_of_interest)

    plot_per_population_provenance_of_interest(
        collated_results, calls, poi, router_provenance_of_interest,
        group_on or group_on_name, current_fig_folder,
        router_pop=router_pop_names)

    plot_population_placement(collated_results, placements,
                              fig_folder=current_fig_folder)

    for folder in run_folders:
        plot_router_provenance(in_folder, folder, router_pop_names,
                               router_provenance_of_interest,
                               current_fig_folder)

    plot_per_population_provenance_of_interest(collated_results, calls, poi,
                                               prov_of_interest,
                                               group_on or group_on_name,
                                               current_fig_folder)

    for folder in run_folders:
        plot_2D_map_for_poi(in_folder, folder,
                            prov_of_interest, router_pop_names,
                            current_fig_folder, placements)


def plot_2D_map_for_poi(folder, selected_sim,
                        provenance_of_interest, router_pop_names, fig_folder, placements):
    write_header("PLOTTING MAPS FOR ALL PRROVENANCE OF INTEREST")
    prov = pd.read_csv(join(join(folder, selected_sim), "structured_provenance.csv"))
    # Filter out router provenance because the logic for plotting those maps is slightly different
    pop_only_prov = prov[~prov['pop'].isin(router_pop_names)]
    filtered_placement = \
        placements[selected_sim]
    try:
        router_provenance = filtered_placement['router_provenance']
    except KeyError:
        traceback.print_exc()
        router_provenance = filtered_placement

    for type_of_provenance in provenance_of_interest:
        #  need to get processor p as well as x y prov_value
        filtered_placement = \
            pop_only_prov[pop_only_prov.prov_name == type_of_provenance][['x', 'y', 'p', 'prov_value']]
        if filtered_placement.shape[0] == 0:
            write_short_msg("NO INFORMATION FOR PROVENANCE", type_of_provenance)
            continue

        # make a new directory for each provenance
        # Check if the results folder exist
        per_prov_dir = os.path.join(fig_folder, type_of_provenance.lower())
        if not os.path.isdir(per_prov_dir) and not os.path.exists(per_prov_dir):
            os.mkdir(per_prov_dir)
        # Plotting bit
        # Fake printing to start things off...
        f = plt.figure(1, figsize=(9, 9), dpi=400)
        plt.close(f)
        # Compute plot order
        plot_order = get_plot_order(router_pop_names)
        plot_display_names = []
        for po in plot_order:
            plot_display_names.append(use_display_name(po))

        magic_constant = 4
        max_x = (router_provenance.x.max() + 1) * magic_constant
        max_y = (router_provenance.y.max() + 1) * magic_constant

        x_ticks = np.arange(0, max_x, magic_constant)[::2]
        x_tick_lables = (x_ticks / magic_constant).astype(int)
        y_ticks = np.arange(0, max_y, magic_constant)[::2]
        y_tick_lables = (y_ticks / magic_constant).astype(int)
        map = np.ones((max_x, max_y)) * np.nan

        for row_index, row in filtered_placement.iterrows():
            x_pos = int(magic_constant * row.x +
                        ((row.p // magic_constant) % magic_constant))
            y_pos = int(magic_constant * row.y +
                        (row.p % magic_constant))
            map[x_pos, y_pos] = row.prov_value

        # crop_point = np.max(np.max(np.argwhere(np.isfinite(map)), axis=0))
        f = plt.figure(1, figsize=(9, 9), dpi=500)
        # plt.matshow(map[:crop_point, :crop_point], interpolation='none')
        im = plt.imshow(map, interpolation='none',
                        cmap=plt.get_cmap('inferno'),
                        extent=[0, max_x, 0, max_y],
                        origin='lower')
        ax = plt.gca()

        plt.xlabel("Chip X coordinate")
        plt.ylabel("Chip Y coordinate")

        plt.xticks(x_ticks, x_tick_lables)
        plt.yticks(y_ticks, y_tick_lables)
        ax.yaxis.set_minor_locator(MultipleLocator(magic_constant))
        ax.xaxis.set_minor_locator(MultipleLocator(magic_constant))

        plt.grid(b=True, which='both', color='k', linestyle='-')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(use_display_name(type_of_provenance))

        save_figure(plt, join(per_prov_dir,
                              "map_of_{}_for_{}".format(type_of_provenance,
                                                        selected_sim)),
                    extensions=['.png', '.pdf'])
        plt.close(f)


def plot_router_provenance(folder, selected_sim, router_pop_names,
                           router_provenance_of_interest, fig_folder):
    write_header("PLOTTING ROUTER INFO AND MAPS")
    prov = pd.read_csv(join(join(folder, selected_sim), "structured_provenance.csv"))
    # Need to filter only info for routers
    # then filter by type of router provenance
    # extract X, Y, prov_value
    # Plot map
    router_only_prov = prov[prov['pop'].isin(router_pop_names)]

    for type_of_provenance in router_provenance_of_interest:
        filtered_placement = \
            router_only_prov[router_only_prov.prov_name == type_of_provenance][['x', 'y', 'prov_value']]
        if filtered_placement.shape[0] == 0:
            write_short_msg("NO INFORMATION FOR PROVENANCE", type_of_provenance)
            continue

        # make a new directory for each provenance
        # Check if the results folder exist
        per_prov_dir = os.path.join(fig_folder, type_of_provenance.lower())
        if not os.path.isdir(per_prov_dir) and not os.path.exists(per_prov_dir):
            os.mkdir(per_prov_dir)
        # Plotting bit
        # Fake printing to start things off...
        f = plt.figure(1, figsize=(9, 9), dpi=400)
        plt.close(f)
        # Compute plot order
        plot_order = get_plot_order(router_pop_names)
        plot_display_names = []
        for po in plot_order:
            plot_display_names.append(use_display_name(po))

        magic_constant = 4
        max_x = (filtered_placement.x.max() + 1) * magic_constant
        max_y = (filtered_placement.y.max() + 1) * magic_constant
        x_ticks = np.arange(0, max_x, magic_constant)[::2]
        # x_tick_lables = np.linspace(0, collated_placements.x.max(), 6).astype(int)
        x_tick_lables = (x_ticks / magic_constant).astype(int)
        y_ticks = np.arange(0, max_y, magic_constant)[::2]
        # y_tick_lables = np.linspace(0, collated_placements.y.max(), 6).astype(int)
        y_tick_lables = (y_ticks / magic_constant).astype(int)
        map = np.ones((max_x, max_y)) * np.nan

        for row_index, row in filtered_placement.iterrows():
            map[
            int(magic_constant * row.x):int(magic_constant * (row.x + 1)),
            int(magic_constant * row.y):int(magic_constant * (row.y + 1))
            ] = row.prov_value

        # crop_point = np.max(np.max(np.argwhere(np.isfinite(map)), axis=0))
        f = plt.figure(1, figsize=(9, 9), dpi=500)
        # plt.matshow(map[:crop_point, :crop_point], interpolation='none')
        im = plt.imshow(map, interpolation='none',
                        cmap=plt.get_cmap('inferno'),
                        extent=[0, max_x, 0, max_y],
                        origin='lower')
        ax = plt.gca()

        plt.xlabel("Chip X coordinate")
        plt.ylabel("Chip Y coordinate")

        plt.xticks(x_ticks, x_tick_lables)
        plt.yticks(y_ticks, y_tick_lables)
        ax.yaxis.set_minor_locator(MultipleLocator(magic_constant))
        ax.xaxis.set_minor_locator(MultipleLocator(magic_constant))

        plt.grid(b=True, which='both', color='k', linestyle='-')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(use_display_name(type_of_provenance))

        save_figure(plt, join(per_prov_dir,
                              "map_of_{}_for_{}".format(type_of_provenance,
                                                        selected_sim)),
                    extensions=['.png', '.pdf'])
        plt.close(f)


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
        # try:
        router_provenance = filtered_placement['router_provenance']
        # except KeyError:
        #     traceback.print_exc()
        #     continue

        placements_per_pop = {x: filtered_placement[x]
                              for x in filtered_placement.keys()
                              if "cell" in x}
        # make a new directory for each provenance
        # Check if the results folder exist
        per_prov_dir = os.path.join(fig_folder, "placements")
        if not os.path.isdir(per_prov_dir) and not os.path.exists(per_prov_dir):
            os.mkdir(per_prov_dir)

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

        magic_constant = 4

        max_x = (router_provenance.x.max() + 1) * magic_constant
        max_y = (router_provenance.y.max() + 1) * magic_constant

        x_ticks = np.arange(0, max_x, magic_constant)[::2]
        x_tick_lables = (x_ticks / magic_constant).astype(int)
        y_ticks = np.arange(0, max_y, magic_constant)[::2]
        y_tick_lables = (y_ticks / magic_constant).astype(int)
        map = np.ones((max_x, max_y)) * np.nan
        for index, pop in enumerate(plot_order):
            curr_pl = placements_per_pop[pop]
            for row_index, row in curr_pl.iterrows():
                x_pos = int(magic_constant * row.x +
                            ((row.p // magic_constant) % magic_constant))
                y_pos = int(magic_constant * row.y +
                            (row.p % magic_constant))
                map[x_pos, y_pos] = index

        uniques = np.unique(map[np.isfinite(map)]).astype(int)
        # crop_point = np.max(np.max(np.argwhere(np.isfinite(map)), axis=0))
        f = plt.figure(1, figsize=(9, 9), dpi=500)
        # plt.matshow(map[:crop_point, :crop_point], interpolation='none')
        im = plt.imshow(map, interpolation='none', vmin=0, vmax=n_plots,
                        cmap=plt.get_cmap('viridis', n_plots),
                        extent=[0, max_x, 0, max_y],
                        origin='lower')
        ax = plt.gca()

        plt.xlabel("Chip X coordinate")
        plt.ylabel("Chip Y coordinate")

        plt.xticks(x_ticks, x_tick_lables)
        plt.yticks(y_ticks, y_tick_lables)
        ax.yaxis.set_minor_locator(MultipleLocator(magic_constant))
        ax.xaxis.set_minor_locator(MultipleLocator(magic_constant))

        plt.grid(b=True, which='both', color='k', linestyle='-')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "5%", pad="3%")
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label("Population")
        cbar.ax.set_yticks(uniques)
        cbar.ax.set_yticklabels(plot_display_names)

        save_figure(plt, join(per_prov_dir,
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


def plot_per_population_provenance_of_interest(
        collated_results, calls, poi, prov,
        group_on, fig_folder, router_pop=None):
    if router_pop is not None:
        write_header("PLOTTING PER ROUTER VALUES FOR PROVENANCE")
    else:
        write_header("PLOTTING PER POPULATION VALUES FOR PROVENANCE")
    sorted_key_list = list(collated_results.keys())
    sorted_key_list.sort()
    f = plt.figure(1, figsize=(9, 9), dpi=400)
    plt.close(f)
    # Look at parameters in calls
    # Check if group_on is a parameter in poi and calls
    # Grab all of the uniques from
    for curr_group in group_on:
        # if curr_g is in POI then we can proceed with plotting things like this
        if curr_group in poi.keys():
            curr_poi = poi[curr_group]
            print("Current parameter of interest is ", curr_group,
                  "with values", curr_poi)
            # sort, just in case
            curr_poi.sort()
            for type_of_prov in prov:
                # Now we need to start collecting values
                # match curr_poi value with collated_results.keys() via calls
                curr_mapping = {val: None for val in curr_poi}
                for curr_poi_val in curr_poi:
                    for call in calls:
                        if call[2][curr_group] == curr_poi_val:
                            # store filename associated with current value
                            match_fname = call[1]
                            if match_fname not in collated_results.keys():
                                continue
                            if type_of_prov not in collated_results[match_fname].keys() and router_pop is not None:
                                filtered_collated_results = {x: None for x in router_pop}
                                for rp in router_pop:
                                    filtered_collated_results[rp] = {
                                        'max': np.nan,
                                        'all': pd.DataFrame([np.nan])}
                            elif type_of_prov in collated_results[match_fname].keys():
                                filtered_collated_results = \
                                    collated_results[match_fname][type_of_prov]
                            else:
                                continue
                            _max_values_per_pop = None
                            if _max_values_per_pop is None:
                                if router_pop is not None:
                                    _max_values_per_pop = \
                                        {x: None for x in filtered_collated_results.keys()
                                         if x in router_pop}
                                else:
                                    _max_values_per_pop = \
                                        {x: None for x in filtered_collated_results.keys()
                                         if "cell" in x}
                            for vp in _max_values_per_pop.keys():
                                _max_values_per_pop[vp] = [list(filtered_collated_results[vp]['all'].values)]

                            # Need to create a list of entries per pop
                            if curr_mapping[curr_poi_val] is None:
                                curr_mapping[curr_poi_val] = _max_values_per_pop
                            else:
                                for k, v in _max_values_per_pop.items():
                                    if v:
                                        curr_mapping[curr_poi_val][k].extend(v)
                # Plotting bit
                plot_order = get_plot_order(_max_values_per_pop.keys())
                n_plots = float(len(plot_order))

                # Convert all 2D arrays of results to numpy arrays
                # Can report intra-trial and inter-trial variability here
                for k in curr_mapping.keys():
                    if curr_mapping[k]:
                        for p in curr_mapping[k].keys():
                            curr_mapping[k][p] = np.array(curr_mapping[k][p])

                f = plt.figure(1, figsize=(9, 9), dpi=400)
                for index, pop in enumerate(plot_order):
                    curr_median = []
                    curr_percentiles = []
                    for k in curr_poi:
                        if curr_mapping[k] and curr_mapping[k][pop].size > 0:
                            merged = np.array(
                                list(itertools.chain.from_iterable(curr_mapping[k][pop])))
                            curr_median.append(np.nanmedian(merged))
                            curr_percentiles.append([np.nanmedian(merged) - np.percentile(merged, 5),
                                                     np.percentile(merged, 95) - np.nanmedian(merged)])
                        else:
                            curr_median.append(np.nan)
                            curr_percentiles.append(np.nan)

                    curr_median = np.array(curr_median)
                    curr_percentiles = np.array(curr_percentiles).T
                    if np.any(np.isfinite(curr_median)):
                        plt.errorbar(curr_poi, curr_median,
                                     yerr=curr_percentiles,
                                     color=color_for_index(index, n_plots),
                                     marker='o',
                                     label=use_display_name(pop),
                                     alpha=0.8)
                        # also print out the values per pop to easily copy and paste
                        write_short_msg(use_display_name(pop), curr_median)
                        write_short_msg(use_display_name(pop), curr_percentiles)
                        # TODO also do some curve fitting for these numbers
                        for deg in [1, 2]:
                            fit_res = polyfit(curr_poi, curr_median, deg)
                            write_short_msg("degree poly {} coeff of "
                                            "determination".format(deg),
                                            fit_res['determination'])
                        fit_res = polyfit(curr_poi, np.log(curr_median), deg)
                        write_short_msg("exp fit coeff of "
                                        "determination",
                                        fit_res['determination'])

                plt.xlabel(use_display_name(curr_group))
                plt.ylabel(use_display_name(type_of_prov))
                ax = plt.gca()
                plt.legend(loc='best')
                plt.tight_layout()
                save_figure(plt, join(fig_folder, "median_{}".format(type_of_prov)),
                            extensions=['.png', '.pdf'])
                if "MAX" in type_of_prov:
                    ax.set_yscale('log')
                    save_figure(plt, join(fig_folder, "log_y_median_{}".format(type_of_prov)),
                                extensions=['.png', '.pdf'])
                plt.close(f)

                f = plt.figure(1, figsize=(9, 9), dpi=400)
                for index, pop in enumerate(plot_order):
                    curr_median = []
                    curr_percentiles = []
                    for k in curr_poi:
                        if curr_mapping[k] and curr_mapping[k][pop].size > 0:
                            merged = np.array(
                                list(itertools.chain.from_iterable(curr_mapping[k][pop])))
                            curr_median.append(np.nanmean(merged))
                            curr_percentiles.append(np.nanstd(merged))
                        else:
                            curr_median.append(np.nan)
                            curr_percentiles.append(np.nan)
                    curr_median = np.array(curr_median)
                    if np.any(np.isfinite(curr_median)):
                        plt.errorbar(curr_poi, curr_median,
                                     yerr=curr_percentiles,
                                     color=color_for_index(index, n_plots),
                                     marker='o',
                                     label=use_display_name(pop),
                                     alpha=0.8)
                        # also print out the values per pop to easily copy and past
                        # also print out the values per pop to easily copy and paste
                        write_short_msg(use_display_name(pop), curr_median)
                        write_short_msg(use_display_name(pop), curr_percentiles)
                        # TODO also do some curve fitting for these numbers
                        for deg in [1, 2]:
                            fit_res = polyfit(curr_poi, curr_median, deg)
                            write_short_msg("degree poly {} coeff of "
                                            "determination".format(deg),
                                            fit_res['determination'])
                        fit_res = polyfit(curr_poi, np.log(curr_median), deg)
                        write_short_msg("exp fit coeff of "
                                        "determination",
                                        fit_res['determination'])

                plt.xlabel(use_display_name(curr_group))
                plt.ylabel(use_display_name(type_of_prov))
                ax = plt.gca()
                plt.legend(loc='best')
                plt.tight_layout()
                save_figure(plt, join(fig_folder, "{}".format(type_of_prov)),
                            extensions=['.png', '.pdf'])

                plt.close(f)
                if "MAX" in type_of_prov:
                    ax.set_yscale('log')
                    save_figure(plt, join(fig_folder, "log_y_median_{}".format(type_of_prov)),
                                extensions=['.png', '.pdf'])
                plt.close(f)


# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

    # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)  # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results


if __name__ == "__main__":
    from spinncer.analysis_argparser import *

    if analysis_args.input and len(analysis_args.input) > 0:
        for in_folder in analysis_args.input:
            # try:
            sweep_provenance_analysis(in_folder, analysis_args.figures_dir,
                                      group_on=analysis_args.group_on,
                                      group_on_name=analysis_args.group_on_name)
            # except Exception as e:
            #     traceback.print_exc()
    else:
        # Constants
        fig_folder = "figures"
        pass
