from spinncer.analysis_common import *
from os.path import join as join


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
        curr_pop_values = filtered_prov[filtered_prov['pop'] == pop].prov_value
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
    # TODO report number of neurons to make sure the networks is correct
    write_short_msg("DETECTED POPULATIONS", pops)

    for type_of_prov in types_of_provenance:
        rep = True if type_of_prov in prov_of_interest else False
        results[type_of_prov] = extract_per_pop_info(prov, type_of_prov, pops,
                                                     report=rep)

    return results, types_of_provenance, prov_of_interest


def sweep_provenance_analysis(in_folder, fig_folder, group_on):
    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)

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
    run_folders = [fc for fc in folder_contents if os.path.isdir(join(in_folder, fc))]
    collated_results = {k: None for k in run_folders}
    types_of_provenance = None
    prov_of_interest = None
    write_short_msg("Number of folders", len(run_folders))
    write_sep()
    for folder in run_folders:
        collated_results[folder], types_of_provenance, \
        prov_of_interest = provenance_csv_analysis(
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


def plot_population_placement():
    pass


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
