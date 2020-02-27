"""
Simulation involving a PyNN script running on either SpiNNaker (through
sPyNNaker) or NEST. The network can either be reduced or full scale, depending
on the connectivity input to it.
"""
import json
# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *

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
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=80)

# Add constraints here
n_neurons_per_core = 64
ss_neurons_per_core = 64
if spinnaker_sim:
    sim.set_number_of_neurons_per_core(sim.IF_cond_exp, n_neurons_per_core)
    sim.set_number_of_neurons_per_core(sim.IF_curr_exp, n_neurons_per_core)
    sim.set_number_of_neurons_per_core(sim.SpikeSourceArray, ss_neurons_per_core)

# Compile stimulus information
stimulus_information = {
    'f_base': args.f_base,
    'f_peak': args.f_peak,
    'stim_times': args.stim_times,
    'stim_radius': args.stim_radius,
    'periodic_stimulus': args.periodic_stimulus,
    'percentage_active_fibers': args.percentage_active_fibers,
}

# Compile json data if it exists
# Load the data from disk
if args.param_json:
    with open(args.param_json, "r") as read_file:
        json_data = json.load(read_file)
else:
    json_data = None

# Instantiate a Cerebellum
cerebellum_circuit = Cerebellum(
    sim, connectivity_filename,
    stimulus_information=stimulus_information,
    reporting=args.no_reports,
    params=json_data,
    skip_projections=args.skip_projections,
    weight_scaling=args.weight_scaling,
    save_conversion_file=args.generate_conversion_constants,
    neuron_model=args.neuron_model
)

# Test various exposed methods
populations = cerebellum_circuit.get_all_populations()

# Set up recordings
cerebellum_circuit.record_all_spikes()
cerebellum_circuit.selectively_record_all(every=n_neurons_per_core)
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
other_recordings = cerebellum_circuit.retrieve_selective_recordings()

# Retrieve final network connectivity
try:
    final_connectivity = cerebellum_circuit.retrieve_final_connectivity()
except:
    # This simulator might not support the way this is done
    final_connectivity = []
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
spike_analysis(results_file=results_file, fig_folder=args.figures_dir)

# Report time taken
print("Results stored in  -- " + filename)

# Report time taken
print("Total time elapsed -- " + str(total_time))
