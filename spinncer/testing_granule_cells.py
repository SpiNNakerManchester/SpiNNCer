"""
Simulation involving a PyNN script running on either SpiNNaker (through
sPyNNaker) or NEST. The network can either be reduced or full scale, depending
on the connectivity input to it.
"""
import json
# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *
from spinncer.utilities.constants import *
import sys

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

# Compile stimulus information
stimulus_information = {
    'f_base': args.f_base,
    'f_peak': args.f_peak,
    'stim_times': [100, 100, 100],
    'stim_radius': np.nan,
    'periodic_stimulus': False,
    'percentage_active_fibers': np.nan,
}

# Compute spikes

n_neurons = 10

# read spikes from file
input_spikes = None
if args.stimulus_from_file is not None:
    # Assembly the dictionary to pass to the cerebellum circuit
    input_spikes = np.load(args.stimulus_from_file, allow_pickle=True)['input_spikes'].ravel()[0]

if input_spikes is None:
    rates = [[0, args.f_peak, 0], ] * n_neurons
    inh_rates = [[0, args.f_base, 0], ] * n_neurons
    starts = [[0, 100, 200], ] * n_neurons
    durations = [[100, 100, 100], ] * n_neurons

    exc_spike_times = create_poisson_spikes(
        n_neurons, rates, starts, durations)

    inh_spike_times = create_poisson_spikes(
        n_neurons, inh_rates, starts, durations)
    input_spikes = {
        'glomerulus': exc_spike_times,
        'golgi': inh_spike_times
    }

# Create Spike Source Arrays
exc_inp = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={
        'spike_times': input_spikes['glomerulus']
    },
    label="glomerulus")

inh_inp = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={
        'spike_times': input_spikes['golgi']
    },
    label="golgi")

# Create LIF population

granule = sim.Population(
    1,
    cellclass=sim.IF_cond_exp,
    cellparams=CELL_PARAMS['granule'],
    label="granule")

# populations = cerebellum_circuit.get_all_populations()
populations = {
    'glomerulus': exc_inp,
    'granule': granule,
    'golgi': inh_inp
}

all_neurons = {
    'glomerulus': n_neurons,
    'granule': 1,
    'golgi': n_neurons
}

if args.generate_conversion_constants:
    np.savez_compressed("conversion_constants",
                        nid_offsets={k: 0 for k in all_neurons.keys()},
                        connectivity=None,
                        all_neurons=all_neurons,
                        connectivity_file=connectivity_filename)
    # expand to also save spikes
    np.savez_compressed(
        "granule_test_stimulus",
        details={
            "argparser": vars(args),
            "git_hash": retrieve_git_commit(),
            "current_time": plt.datetime.datetime.now().strftime("%H:%M:%S_%d/%m/%Y"),
            "n_neurons_per_core": 255,
            "ss_neurons_per_core": 255,
        },
        input_spikes=input_spikes,
        **input_spikes)
    print("Saved conversion_constants.npz. Exiting...")
    sys.exit(0)

# Create Projections
glom_grc = sim.Projection(
    exc_inp,  # pre-synaptic population
    granule,  # post-synaptic population
    # connector includes (source, target, weight, delay)
    sim.AllToAllConnector(),
    synapse_type=sim.StaticSynapse(
        weight=CONNECTIVITY_MAP['glom_grc']['weight'],
        delay=args.timestep),
    receptor_type="excitatory",  # inh or exc
    label="glom_grc")  # label for connection

goc_grc = sim.Projection(
    exc_inp,  # pre-synaptic population
    granule,  # post-synaptic population
    # connector includes (source, target, weight, delay)
    sim.AllToAllConnector(),
    synapse_type=sim.StaticSynapse(
        weight=CONNECTIVITY_MAP['goc_grc']['weight'],
        delay=args.timestep*6),
    receptor_type="inhibitory",  # inh or exc
    label="goc_grc")  # label for connection

projections = {
    'glom_grc': glom_grc,
    'goc_grc': goc_grc
}

# Set up recordings
# cerebellum_circuit.record_all_spikes()
# cerebellum_circuit.selectively_record_all(every=1)
for label, pop in populations.items():
    print("Enabling spike recordings for ", label, "...")
    pop.record(['spikes'])
    if label == "granule":
        print("Selectively recording gsyn and v for ", label)
        ps = pop.size
        pop.record(['gsyn_inh', 'gsyn_exc', 'v'])
    else:
        print("Skipping selective recording for", label, "...")

recorded_spikes = {}
other_recordings = {}

# Record simulation start time (wall clock)
sim_start_time = plt.datetime.datetime.now()
current_error = None
# Run the simulation
try:
    sim.run(args.simtime)  # ms
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
    if label in ["glomerulus", "golgi"]:
        print("Skipping selective recording for", label, "...")
        continue
    print("Retrieving recordings for ", label, "...")
    other_recordings[label] = {}
    other_recordings[label]['gsyn_inh'] = pop.get_data(['gsyn_inh'])
    other_recordings[label]['gsyn_exc'] = pop.get_data(['gsyn_exc'])
    other_recordings[label]['v'] = pop.get_data(['v'])

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
    filename = "granule_cell_test" + str(suffix)

if current_error:
    filename = "error_" + filename

# Check if the results folder exist
if not os.path.isdir(args.result_dir) and not os.path.exists(args.result_dir):
    os.mkdir(args.result_dir)

# Retrieve simulation parameters for provenance tracking and debugging purposes
sim_params = {
    "argparser": vars(args),
    "git_hash": retrieve_git_commit(),
    "run_end_time": end_time.strftime("%H:%M:%S_%d/%m/%Y"),
    "wall_clock_script_run_time": str(total_time),
    "wall_clock_sim_run_time": str(sim_total_time),
    "n_neurons_per_core": 1,
    "ss_neurons_per_core": 1,
}

# Save results to file in [by default] the `results/' directory
results_file = os.path.join(args.result_dir, filename)
np.savez_compressed(results_file,
                    simulation_parameters=sim_params,
                    all_spikes=recorded_spikes,
                    other_recordings=other_recordings,
                    all_neurons=all_neurons,
                    final_connectivity=final_connectivity,
                    initial_connectivity=initial_connectivity,
                    stimulus_params=stimulus_information,
                    simtime=args.simtime,
                    json_data=None,
                    conn_params=CONNECTIVITY_MAP,
                    cell_params=CELL_PARAMS,
                    )

# Appropriately end the simulation
sim.end()

# Analysis time!
spike_analysis(results_file=results_file, fig_folder=args.figures_dir,
               worst_case=False, delay_sensitive=False)

# Report time taken
print("Results stored in  -- " + filename)

# Report time taken
print("Total time elapsed -- " + str(total_time))
