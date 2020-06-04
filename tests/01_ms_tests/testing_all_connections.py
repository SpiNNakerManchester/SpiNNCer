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
    try:
        # this might be deprecated soon
        import spynnaker8 as sim
    except ImportError:
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

# Record SCRIPT start time (wall clock)
start_time = plt.datetime.datetime.now()

connectivity_filename = args.dataset or DEFAULT_DATASET
connectivity_filename = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'datasets', connectivity_filename)
# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=1,
          timescale=args.timescale)


RB_LEFT_SHIFTS = {
    'golgi': [0, 0],
    'granule': [0, 0],
    'purkinje': [-3, 0],
    'basket': [-2, 0],
    'stellate': [-2, 0],
    'dcn': [-5, 0]
}

n_neurons = 10

# Compute spikes
input_spike = 10 - args.timestep

# Create Spike Source Arrays
single_spike_source = sim.Population(
    1,
    cellclass=sim.SpikeSourceArray,
    cellparams={
        'spike_times': [input_spike]
    },
    label="glomerulus")


# Create a LIF population per projection
populations = {}
projections = {}
additional_parameters = {}

initial_connectivity = {}
for conn_name, conn_params in CONNECTIVITY_MAP.items():
    cell_params = CELL_PARAMS[conn_params['post']]
    if str.lower(args.simulator) in ["spinnaker", "spynnaker"]:
        additional_parameters = {
            "additional_parameters":{
                "rb_left_shifts": RB_LEFT_SHIFTS[conn_params['post']],
                "n_steps_per_timestep": args.loops_grc
            }
        }
    print("ADDITIONAL PARAMETERS:", additional_parameters)
    # Create population
    if args.disable_i_offset and "i_offset" in cell_params.keys():
        del cell_params['i_offset']
        print("Removed i_offset from ", conn_params['post'])
    cell_to_test = sim.Population(
        1,
        cellclass=sim.IF_cond_exp,
        cellparams=cell_params,
        label=conn_name,
        **additional_parameters)
    # initialise V
    cell_to_test.initialize(v=cell_params['v_rest'])
    populations[conn_name] = cell_to_test
    # Create projection
    conn_to_test = sim.Projection(single_spike_source, cell_to_test,
                                  sim.OneToOneConnector(),
                                  synapse_type=sim.StaticSynapse(
                                      weight=conn_params['weight'],
                                      delay=args.timestep),
                                  receptor_type="excitatory" if conn_params['weight'] > 0 else "inhibitory",
                                  label=conn_name)
    projections[conn_name] = conn_to_test
    # Enable recordings
    cell_to_test.record(['spikes', 'gsyn_inh', 'gsyn_exc', 'v'])

all_neurons = {k:1 for k in populations.keys()}
recorded_spikes = {}
other_recordings = {}

# Record simulation start time (wall clock)
sim_start_time = plt.datetime.datetime.now()
current_error = None
# Run the simulation
try:
    sim.run(100)  # ms
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
        print("Retrieving recordings for ", label, "...")
        recorded_spikes[label] = pop.get_data(['spikes'])
# other_recordings = cerebellum_circuit.retrieve_selective_recordings()
other_recordings = {}
for label, pop in populations.items():
    print("Retrieving recordings for ", label, "...")
    other_recordings[label] = {}
    other_recordings[label]['gsyn_inh'] = pop.get_data(['gsyn_inh']).filter(name='gsyn_inh')[0].magnitude.ravel()
    other_recordings[label]['gsyn_exc'] = pop.get_data(['gsyn_exc']).filter(name='gsyn_exc')[0].magnitude.ravel()
    other_recordings[label]['v'] = pop.get_data(['v']).filter(name='v')[0].magnitude.ravel()

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
excel_filename = "{}_{}_loops_{}_testing_all_connections.xlsx".format(
    args.simulator,
    args.loops_grc,
    "WITH_ioffset" if not args.disable_i_offset else "WITHOUT_ioffset"
)

writer = pd.ExcelWriter(
    os.path.join(args.result_dir, excel_filename),
    engine='xlsxwriter')

ordered_projections = list(other_recordings.keys())
ordered_projections.sort()
for key in ordered_projections:
    value = other_recordings[key]
    df = pd.DataFrame.from_dict(value)
    df = df[['v', 'gsyn_exc', 'gsyn_inh']]
    df.to_excel(writer, sheet_name=key)

writer.save()

# Report time taken
print("Total time elapsed -- " + str(total_time))
