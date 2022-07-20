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
from spinncer.utilities import create_poisson_spikes

spinnaker_sim = False
if str.lower(args.simulator) in ["spinnaker", "spynnaker"]:
    import pyNN.spynnaker as sim
    spinnaker_sim = True
elif str.lower(args.simulator) in ["nest"]:
    import pyNN.nest as sim
else:
    raise ValueError("Simulator " + str.lower(args.simulator) +
                     "unrecognised!")
# provenance utility
from spinncer.utilities.provenance import retrieve_git_commit
# analysis functionsc
from spinncer.cerebellum_analysis import *
import pylab as plt
import os
import pandas as pd
import xlsxwriter
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
          timescale=args.timescale)

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
if str.lower(args.simulator) in ["spinnaker", "spynnaker"]:
    # additional_params = {"rb_left_shifts": [0, 0]}
    additional_parameters = {
        "additional_parameters": {
            "rb_left_shifts": [0, 0],
            "n_steps_per_timestep": args.loops_grc
        }
    }
    print("ADDITIONAL PARAMETERS:", additional_parameters)

for conn_name, conn_params in CONNECTIVITY_MAP.items():
    cell_params = CELL_PARAMS[conn_params['post']]
    # Create population
    if args.disable_i_offset and "i_offset" in cell_params.keys():
        del cell_params['i_offset']
        print("Removed i_offset from ", conn_params['post'])
    # flag if connections is exc
    is_conn_exc = conn_params['weight'] > 0
    is_projection_exc[conn_name] = 0 if is_conn_exc else 1
    cell_to_test = sim.Population(
        1,
        cellclass=sim.IF_cond_exp,
        cellparams=cell_params,
        label=args.population,
        **additional_parameters)
    # initialise V
    cell_to_test.initialize(v=cell_params['v_rest'])
    populations[conn_name] = cell_to_test
    curr_weight = np.abs(conn_params['weight'])  # NEST is awful.
    if args.test_max_weight:
        curr_weight *= EXPECTED_MAX_SPIKES[conn_name]
    # Create projection
    conn_to_test = sim.Projection(all_spike_sources[conn_name], cell_to_test,
                                  sim.OneToOneConnector(),
                                  synapse_type=sim.StaticSynapse(
                                      weight=curr_weight,
                                      delay=args.timestep),
                                  receptor_type="excitatory" if conn_params['weight'] > 0 else "inhibitory",
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
        recorded_spikes[label] = pop.get_data(['spikes']).segments[0].spiketrains[0]
        print("Spikes: ", recorded_spikes[label])
# other_recordings = cerebellum_circuit.retrieve_selective_recordings()
other_recordings = {}
gsyn_labels = ['gsyn_inh', 'gsyn_exc']
for label, pop in populations.items():
    print("Retrieving recordings for ", label, "...")
    other_recordings[label] = {}
    curr_gsyn_type = gsyn_labels[is_projection_exc[label]]

    other_recordings[label]['gsyn'] = np.array(
        pop.get_data([curr_gsyn_type]).filter(name=curr_gsyn_type))[0].T.ravel()
    other_recordings[label]['v'] = np.array(
        pop.get_data(['v']).segments[0].filter(name='v'))[0].T.ravel()
    other_recordings[label]['current'] = np.ones(other_recordings[label]['v'].shape) * np.nan

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
initial_connectivity = []

# Save results
suffix = end_time.strftime("_%H%M%S_%d%m%Y")
if args.filename:
    filename = args.filename
else:
    filename = args.population + "_cell_test" + str(suffix)

if current_error:
    filename = "error_" + filename

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
    spike_hit_vector = np.zeros(int((simtime) * 10)+1)
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
# Report time taken
print("Results stored in  -- " + filename)

# Report time taken
print("Total time elapsed -- " + str(total_time))
