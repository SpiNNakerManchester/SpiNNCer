"""
Simulation involving a PyNN script running on either SpiNNaker (through
sPyNNaker) or NEST. The network consists only of populations with single
neurons and no connections. This script is meant for testing the behaviour
of single cells using paramters to be used in the full-scale
cerebellum_experiment.py
"""
from spinncer.spinncer_argparser import *
import sys
import json

# import sPyNNaker
# import simulator
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
from spinncer.cerebellum import Cerebellum
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
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=6.4)

# Add constraints here
n_neurons_per_core = 1
ss_neurons_per_core = 1
if spinnaker_sim:
    sim.set_number_of_neurons_per_core(sim.IF_cond_exp, n_neurons_per_core)
    sim.set_number_of_neurons_per_core(sim.IF_curr_exp, n_neurons_per_core)
    sim.set_number_of_neurons_per_core(sim.SpikeSourceArray, ss_neurons_per_core)

# Compile stimulus information

# Modify data here because of the restricted meaning of the experiment
args.stim_times = [100.0, 100.0, 300.0]
args.simtime = 500
args.skip_projections = True

stimulus_information = {
    'f_base': args.f_base,
    'f_peak': args.f_peak,
    'stim_times': args.stim_times,
    'stim_radius': args.stim_radius,
    'periodic_stimulus': args.periodic_stimulus,
}

# Compile json data if it exists
# Load the data from disk
if args.param_json:
    with open(args.param_json, "r") as read_file:
        json_data = json.load(read_file)
else:
    json_data = None

# read spikes from file
input_spikes = None
if args.stimulus_from_file is not None:
    # Assembly the dictionary to pass to the cerebellum circuit
    input_spikes = np.load(args.stimulus_from_file, allow_pickle=True)['input_spikes'].ravel()[0]

# Instantiate a Cerebellum with *1 neuron per population without projections.*
# We add those outside here.
cerebellum_circuit = Cerebellum(
    sim, connectivity_filename,
    stimulus_information=stimulus_information,
    reporting=args.no_reports,
    params=json_data,
    force_number_of_neurons=1,
    skip_projections=True,
    weight_scaling=args.weight_scaling,
    neuron_model=args.neuron_model,
    input_spikes=input_spikes,
    rb_left_shifts=args.rb_left_shifts
)

if args.generate_conversion_constants:
    np.savez_compressed("conversion_constants",
                        nid_offsets=cerebellum_circuit.nid_offset,
                        connectivity=cerebellum_circuit.connections,
                        all_neurons=cerebellum_circuit.number_of_neurons,
                        connectivity_file=connectivity_filename)
    # expand to also save spikes
    np.savez_compressed(
        "stimulus",
        details={
            "argparser": vars(args),
            "git_hash": retrieve_git_commit(),
            "current_time": plt.datetime.datetime.now().strftime("%H:%M:%S_%d/%m/%Y"),
            "n_neurons_per_core": n_neurons_per_core,
            "ss_neurons_per_core": ss_neurons_per_core,
        },
        input_spikes=cerebellum_circuit.stimulus,
        **cerebellum_circuit.stimulus)
    print("Saved conversion_constants.npz. Exiting...")
    sys.exit(0)
# Assemble the stimulus
# Copy-paste from cell.pptz
stim_array = [100., 110., 120., 130., 140., 150., 160., 170., 180., 190.]
# Create the stimulus spike source array
stimulus_pop = sim.Population(
    1,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': stim_array},
    label="Single cell stimulus")

# Retrieve all of the populations in the cerebellum
populations = cerebellum_circuit.get_all_populations()

# Form 1:1 connections to each of them
# The following values for weight and delay are provided by BM 18.02.2020
WEIGHT = 0.01
DELAY = 1.0
stim_to_pop_proj = []
for pop_name, pop in populations.items():
    if pop and 'glom' not in pop_name:
        stim_to_pop_proj.append(sim.Projection(
            stimulus_pop, pop, sim.OneToOneConnector(),
            sim.StaticSynapse(weight=WEIGHT, delay=DELAY),
            receptor_type="excitatory",
            label="conn_from_stimulus_to_" + pop_name)
        )

# Set up recordings
cerebellum_circuit.record_all_spikes()
cerebellum_circuit.selectively_record_all(every=n_neurons_per_core)
# Record spikes from the SpikeSourceArray
stimulus_pop.record(['spikes'])
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
recorded_spikes = cerebellum_circuit.retrieve_all_recorded_spikes()
recorded_spikes['glomerulus'] = stimulus_pop.get_data(['spikes'])
other_recordings = cerebellum_circuit.retrieve_selective_recordings()

# Retrieve final network connectivity
try:
    final_connectivity = cerebellum_circuit.retrieve_final_connectivity()
except:
    # This simulator might not support the way this is done
    final_connectivity = []
    traceback.print_exc()
initial_connectivity = cerebellum_circuit.connections

# Save results
suffix = end_time.strftime("_%H%M%S_%d%m%Y")
if args.filename:
    filename = args.filename
else:
    filename = "spinncer_experiment" + str(suffix)

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
    "n_neurons_per_core": n_neurons_per_core,
    "ss_neurons_per_core": ss_neurons_per_core,
}

# Save results to file in [by default] the `results/' directory
results_file = os.path.join(args.result_dir, filename)
np.savez_compressed(results_file,
                    simulation_parameters=sim_params,
                    all_spikes=recorded_spikes,
                    other_recordings=other_recordings,
                    all_neurons=cerebellum_circuit.number_of_neurons,
                    final_connectivity=final_connectivity,
                    initial_connectivity=initial_connectivity,
                    stimulus_params=stimulus_information,
                    simtime=args.simtime,
                    json_data=json_data,
                    conn_params=cerebellum_circuit.retrieve_conn_params(),
                    cell_params=cerebellum_circuit.retrieve_cell_params(),
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
