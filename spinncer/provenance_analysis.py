from spinncer.analysis_common import *
from os.path import join as join


def provenance_csv_analysis(in_folder, fig_folder):
    write_header("Reading provenances in folder " + in_folder)
    prov = pd.read_csv(join(in_folder, "structured_provenance.csv"))
    pops = prov['pop'].unique()
    types_of_provenance = prov['prov_name'].unique()
    max_spikes_in_a_tick = {}
    max_dmas_in_a_tick = {}

    # TODO report number of neurons to make sure the networks is correct
    write_short_msg("DETECTED POPULATIONS", pops)

    # check that no cores have a positive value for 'max_number_of_times_timer_tic_over_ran'
    if 'max_number_of_times_timer_tic_over_ran' in types_of_provenance:
        overran = prov['prov_name'] == 'max_number_of_times_timer_tic_over_ran'
        where_overran = prov[overran].prov_value > 0
        if any(where_overran):
            print(Fore.RED + "SOME CORES OVERRAN THE TIMER TICK")
            print("\t" + prov[where_overran].prov_name.unique(),
                  Style.RESET_ALL)
    else:
        print(Fore.GREEN + "max_number_of_times_timer_tic_over_ran not in file",
              Style.RESET_ALL)

    # TODO report 'MAX_SPIKES_IN_A_TICK' for each population (avg and max)
    if 'MAX_SPIKES_IN_A_TICK' in types_of_provenance:
        max_tick_filter = prov['prov_name'] == 'MAX_SPIKES_IN_A_TICK'
        max_tick_prov = prov[max_tick_filter]
        print("MAX SPIKE IN A TICK FOR POPULATIONS:")
        for pop in pops:
            curr_pop_values = max_tick_prov[max_tick_prov['pop'] == pop].prov_value
            print("\t{:25} - avg {:10.2f} max {:10.2f}".format(
                pop, curr_pop_values.mean(), curr_pop_values.max()))
    else:
        print(Fore.GREEN + "MAX_SPIKES_IN_A_TICK not in file",
              Style.RESET_ALL)

    # TODO report 'MAX_DMAS_IN_A_TICK' for each population (avg and max)


    # TODO report 'Times_synaptic_weights_have_saturated' for each population (avg and max)

    return max_spikes_in_a_tick


def sweep_provenance_analysis(in_folder, fig_folder):
    # Check if the folders exist
    if not os.path.isdir(fig_folder) and not os.path.exists(fig_folder):
        os.mkdir(fig_folder)

    write_header("Analysing provenances in folder " + in_folder)
    folder_contents = os.listdir(in_folder)
    batch_infos = [fc for fc in folder_contents if "batch_" in fc]
    if len(batch_infos) > 1:
        raise FileExistsError("There are multiple batch info .npz archives "
                              "in the current directory!")
    batch_info = np.load(join(in_folder, batch_infos[0]))
    run_id = "@" + batch_infos[0].split("_")[1].split(".")[0]
    if batch_info:
        write_header("Batch info contained in " + batch_infos[0])
        write_short_msg("FILES", batch_info.files)
        poi = batch_info['parameters_of_interest'].ravel()[0]
        write_short_msg("PARAMETERS", poi)
    run_folders = [fc for fc in folder_contents if run_id in fc]
    collated_results = {k: None for k in run_folders}
    write_short_msg("Number of folders", len(run_folders))
    write_sep()
    for folder in run_folders:
        collated_results[folder] = provenance_csv_analysis(
            join(in_folder, folder), fig_folder)


if __name__ == "__main__":
    from spinncer.analysis_argparser import *

    if analysis_args.input and len(analysis_args.input) > 0:
        for in_folder in analysis_args.input:
            sweep_provenance_analysis(in_folder, analysis_args.figures_dir)
    else:
        # Constants
        fig_folder = "figures"
        pass
