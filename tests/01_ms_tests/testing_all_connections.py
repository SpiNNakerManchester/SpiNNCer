"""
Simulation involving a PyNN script running on either SpiNNaker (through
sPyNNaker) or NEST. Test all individual connections in the network and their
effect on each relevant post-synaptic neuron
"""
import json
# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *
from spinncer.utilities.constants import *
import sys
from spinncer.utilities.utils import floor_spike_time

# import sPyNNaker
# import simulator
from spinncer.utilities import create_poisson_spikes, round_to_nearest_accum

spinnaker_sim = False
if str.lower(args.simulator) in ["spinnaker", "spynnaker"]:
    import pyNN.spynnaker as sim
    spinnaker_sim = True
elif str.lower(args.simulator) in ["nest"]:
    import pyNN.nest as sim
else:
    raise ValueError("Simulator " + str.lower(args.simulator) +
                     "unrecognised!")
# analysis functionsc
from spinncer.cerebellum_analysis import *
import pylab as plt
import os
import pandas as pd
import numpy as np
import traceback

# Fail if passing in incorrect parameters
if args.test_max_weight and args.test_max_spikes:
    raise AttributeError("Can either test maximum weight effect on neuron (single spike, weight scaled up) "
                         "or maximum number of spikes (number of spikes scaled up, weight as original).")

# Record SCRIPT start time (wall clock)
start_time = plt.datetime.datetime.now()

connectivity_filename = args.dataset or DEFAULT_DATASET
connectivity_filename = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'datasets', connectivity_filename)
# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=1,
          timescale=args.timescale,
          spike_precision=args.nest_grid  # NEST Spike precision
          )

simtime = args.simtime

EXPECTED_MAX_SPIKES = {
    # Empirical values observed in simulation with a 200 Hz input with stim_radius = 130
    'aa_goc': 31,
    'aa_pc': 24,
    'bc_pc': 7,
    'gj_bc': 4,
    'gj_goc': 25,
    'gj_sc': 6,
    'glom_dcn': 6,
    'glom_goc': 8,
    'glom_grc': 5,
    'goc_grc': 4,
    'pc_dcn': 6,
    'pf_bc': 64,
    'pf_goc': 64,
    'pf_pc': 729,
    'pf_sc': 60,
    'sc_pc': 9,
}

# Number of pairs
n_pops = len(list(CONNECTIVITY_MAP.keys()))

# Compute spikes
input_spike_time = 10 - args.timestep

# Create Spike Source Arrays
if not args.test_max_spikes:
    single_spike_source = sim.Population(
        2,
        cellclass=sim.SpikeSourceArray,
        cellparams={
            'spike_times': [input_spike_time]
        },
        label="Spike source")

all_spike_sources = {}
for conn_label in CONNECTIVITY_MAP.keys():
    if not args.test_max_spikes:
        all_spike_sources[conn_label] = single_spike_source
    else:
        all_spike_sources[conn_label] = sim.Population(
            2,
            cellclass=sim.SpikeSourceArray,
            cellparams={
                'spike_times': [input_spike_time] * EXPECTED_MAX_SPIKES[conn_label]
            },
            label="Spike source for {}".format(conn_label))

is_projection_exc = {}

# Create a LIF population per projection
populations = {}
projections = {}
additional_parameters = {}
per_pop_r_mem = {}

canonical_rbls = RMEM_RBLS if args.r_mem else VANILLA_RBLS
initial_connectivity = {}
for conn_name, conn_params in CONNECTIVITY_MAP.items():
    print("-" * 80)
    print("CONNECTION ", conn_name)
    cell_params = copy.deepcopy(CELL_PARAMS[conn_params['post']])

    curr_weight = conn_params['weight']

    # pre-multiply membrane Resistance into weight and i_offset?
    if args.r_mem:
        r_mem = cell_params['tau_m'] / cell_params['cm']
        print("R_mem", r_mem)
        # adjust i_offset
        if 'i_offset' in cell_params.keys():
            print("original i_offset =", cell_params['i_offset'])
            cell_params['i_offset'] *= r_mem
            print("r_mem * i_offset =", cell_params['i_offset'])
        # adjust weight
        curr_weight *= r_mem
        per_pop_r_mem[conn_name] = r_mem
    else:
        per_pop_r_mem[conn_name] = 1.0

    if args.test_max_weight:
        curr_weight *= EXPECTED_MAX_SPIKES[conn_name]
    # flag if connections is exc
    is_conn_exc = curr_weight > 0
    if is_conn_exc:
        weight_m = [np.abs(curr_weight), 0]
    else:
        weight_m = [0, np.abs(curr_weight)]
    is_projection_exc[conn_name] = 0 if is_conn_exc else 1

    if args.real_rbls:
        # Retrieve the RB LS for the post-synaptic population
        rb_ls = canonical_rbls[conn_params['post']]
        print("Using cannonical RB LS for population", args.population, "with values", rb_ls)
    else:
        if args.r_mem:
            if not args.test_max_spikes:
                rb_ls = np.asarray([np.ceil(np.log2(
                    np.abs(curr_weight))), ] * 2)
            else:
                rb_ls = np.asarray([np.ceil(np.log2(
                    np.abs(curr_weight * EXPECTED_MAX_SPIKES[conn_name]))), ] * 2)

        else:
            if not args.test_max_spikes:
                rb_ls = np.asarray([np.ceil(np.log2(
                    np.abs(curr_weight * (2 ** 5)))), ] * 2)
            else:
                rb_ls = np.asarray([np.ceil(np.log2(
                    np.abs(curr_weight * (2 ** 5) * EXPECTED_MAX_SPIKES[conn_name]))), ] * 2)
        print("Compute RB_LS:", rb_ls)
        rb_ls = np.clip(rb_ls, 0, 16)
        print("Clipped RB_LS:", rb_ls)
    additional_parameters = {
        "additional_parameters": {
            "rb_left_shifts": rb_ls,
            "n_steps_per_timestep": args.loops_grc
        }
    }
    print("ADDITIONAL PARAMETERS:", conn_params['post'], additional_parameters)

    # Remove i_offset if we don't want it in this experiment
    if args.disable_i_offset and "i_offset" in cell_params.keys():
        del cell_params['i_offset']
        print("Removed i_offset from ", conn_params['post'])
    # Create population
    cell_to_test = sim.Population(
        2,
        cellclass=sim.IF_cond_exp,
        cellparams=cell_params,
        label=conn_name,
        **additional_parameters)
    # initialise V
    cell_to_test.initialize(v=cell_params['v_rest'])
    populations[conn_name] = cell_to_test
    print(conn_name, "will have weight set to", weight_m)
    # Create projection
    conn_to_test = sim.Projection(all_spike_sources[conn_name], cell_to_test,
                                  sim.OneToOneConnector(),
                                  synapse_type=sim.StaticSynapse(
                                      weight=weight_m,
                                      delay=args.timestep),
                                  receptor_type="excitatory" if is_conn_exc else "inhibitory",
                                  label=conn_name)
    projections[conn_name] = conn_to_test
    # Enable recordings
    cell_to_test.record(['spikes', 'gsyn_inh', 'gsyn_exc', 'v'])

all_neurons = {k: 1 for k in populations.keys()}
recorded_spikes = {}
other_recordings = {}

# Record simulation start time (wall clock)
sim_start_time = plt.datetime.datetime.now()
current_error = None
# Run the simulation
try:
    sim.run(simtime)  # ms
except Exception as e:
    print("An exception occurred during execution!")
    traceback.print_exc()
    current_error = e

# Compute time taken to reach this point
end_time = plt.datetime.datetime.now()
total_time = end_time - start_time
sim_total_time = end_time - sim_start_time

# Retrieve recordings
# recorded_spikes = cerebellum_circuit.retrieve_all_recorded_spikes()
recorded_spikes = {}
for label, pop in populations.items():
    if pop is not None:
        print("Retrieving spikes for ", label, "...")
        recorded_spikes[label] = pop.get_data(['spikes']).segments[0].spiketrains[is_projection_exc[label]]
        print("Spikes: ", recorded_spikes[label])
# other_recordings = cerebellum_circuit.retrieve_selective_recordings()
other_recordings = {}
for label, pop in populations.items():
    print("Retrieving recordings for ", label, "...")
    other_recordings[label] = {}

    other_recordings[label]['current'] = np.array(
        pop.get_data(['gsyn_inh']).filter(name='gsyn_inh'))[0].T[is_projection_exc[label]] / per_pop_r_mem[label]

    other_recordings[label]['gsyn'] = np.array(
        pop.get_data(['gsyn_exc']).filter(name='gsyn_exc'))[0].T[is_projection_exc[label]].ravel() / per_pop_r_mem[
                                          label]

    other_recordings[label]['v'] = np.array(
        pop.get_data(['v']).segments[0].filter(name='v'))[0].T[is_projection_exc[label]]

# Retrieve final network connectivity
try:
    # final_connectivity = cerebellum_circuit.retrieve_final_connectivity()
    final_connectivity = {}
    for label, p in projections.items():
        if p is None:
            print("Projection", label, "is not implemented!")
            continue
        print("Retrieving connectivity for projection ", label, "...")
        try:
            final_connectivity[label] = \
                np.array(p.get(('weight', 'delay'),
                               format="list")._get_data_items())
        except Exception as e:
            print("Careful! Something happened when retrieving the "
                  "connectivity:", e, "\nRetrying...")
            final_connectivity[label] = \
                np.array(p.get(('weight', 'delay'), format="list"))
except:
    # This simulator might not support the way this is done
    final_connectivity = []
    traceback.print_exc()

print("MAX weight per projection")
print("-" * 80)
conn_dict = {}
print("{:27} -> SpiNNaker weight \t|\t Prescribed weight".format("Conn name"))
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
    mean = np.max(conn[:, 2]) / per_pop_r_mem[key]
    if args.test_max_weight:
        mean = mean / EXPECTED_MAX_SPIKES[key]
    # replace with percentage of difference
    original_conn = np.abs(CONNECTIVITY_MAP[key]["weight"])
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
            CONNECTIVITY_MAP[key]["weight"], proportion))

# Check if the results folder exist
if not os.path.isdir(args.result_dir) and not os.path.exists(args.result_dir):
    os.mkdir(args.result_dir)

# Appropriately end the simulation
sim.end()

# save required csv
excel_filename = "testing_all_connections_{}_{}_loops".format(
    args.simulator,
    args.loops_grc,
)
if args.r_mem:
    excel_filename += "_" + "r_mem"
if args.test_max_spikes:
    excel_filename += "_" + "max_spikes"
if args.test_max_weight:
    excel_filename += "_" + "single_max_weight"
if args.disable_i_offset:
    excel_filename += "_" + "WITHOUT_IOFFSET"

if args.suffix:
    excel_filename += "_" + args.suffix
excel_filename += ".xlsx"

writer = pd.ExcelWriter(
    os.path.join(args.result_dir, excel_filename),
    engine='xlsxwriter')

ordered_projections = list(other_recordings.keys())
ordered_projections.sort()
for key in ordered_projections:
    value = other_recordings[key]

    # Create a bool vector with an entry per time step. 0 = no spike, 1 = spike
    spike_hit_vector = np.zeros(int(simtime * 10))
    curr_spikes = np.asarray(recorded_spikes[key]).ravel()
    curr_timestep = (curr_spikes * 10).astype(int)
    spike_hit_vector[curr_timestep] = 1

    # add spikes to value dict
    value["spikes"] = spike_hit_vector

    df = pd.DataFrame.from_dict(value)
    # just order columns alphabetically
    ordered_variables = list(value.keys())
    ordered_variables.sort()
    df = df[ordered_variables]
    df.to_excel(writer, sheet_name=key)

writer.save()

pp(per_pop_r_mem)

# Report time taken
print("Total time elapsed -- " + str(total_time))
