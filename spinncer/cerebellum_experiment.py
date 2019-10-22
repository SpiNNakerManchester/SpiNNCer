# import sPyNNaker
from spinncer.cerebellum import Cerebellum
import numpy as np

try:
    # this might be deprecated soon
    import spynnaker8 as sim
except:
    import pyNN.spynnaker as sim
# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *
# provenance utility
from spinncer.utilities.provenance import retrieve_git_commit
import pylab as plt
import os

# Record SCRIPT start time (wall clock)
start_time = plt.datetime.datetime.now()
connectivity_filename = 'datasets/scaffold_detailed__158.0x158.0_v3.hdf5'

# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=10)

cerebellum_circuit = Cerebellum(sim, connectivity_filename,
                                reporting=args.no_reports,
                                skip_projections=args.skip_projections)

# Test various exposed methods
populations = cerebellum_circuit.get_all_populations()
assert (len(populations) == 7)
input_populations = cerebellum_circuit.get_circuit_inputs()
assert (len(input_populations) == 1)
output_populations = cerebellum_circuit.get_circuit_outputs()
assert (len(output_populations) == 1)

# Add stimulus to the network -- GLOM
# number of neurons for input is number of neurons in GLOM
n_inputs = cerebellum_circuit.number_of_neurons['glomerulus']
# convert stimulation times to numpy array
stim_times = np.asarray(args.stim_times)
# compute number of rate changes required
number_of_slots = int(stim_times.size)
# compute the time at which individual rates take effect
stim_start = np.cumsum(np.append([0], stim_times))[:number_of_slots]
starts = np.ones((n_inputs, number_of_slots)) * stim_start
# compute the duration (in ms) for which individual rates are active
durations = np.ones((n_inputs, number_of_slots)) * stim_times
# compute the individual rates (in Hz) for each slot
rates = np.ones((n_inputs, number_of_slots)) * \
        np.asarray([args.f_base, args.f_base + args.f_peak, args.f_base])
# Add a variable-rate poisson spike source to the network
stimulus_params = {
    'rates': rates,
    'starts': starts,
    'durations': durations,
}
stimulus = sim.Population(n_inputs, sim.extra_models.SpikeSourcePoissonVariable,
                          stimulus_params, label='stimulus population')
# Connect stimulus to the relevant populations (here, granular)
# The connection has the same parameters as glom_grc
sim.Projection(stimulus, cerebellum_circuit.granule,
               sim.FromListConnector(cerebellum_circuit.connections['glom_grc']),
               label="stimulus projection")


# Set up recordings
cerebellum_circuit.record_all_spikes()

# Record simulation start time (wall clock)
sim_start_time = plt.datetime.datetime.now()

# Run the simulation
sim.run(args.simtime)  # ms

# Compute time taken to reach this point
end_time = plt.datetime.datetime.now()
total_time = end_time - start_time
sim_total_time = end_time - sim_start_time

# Retrieve recordings
recorded_spikes = cerebellum_circuit.retrieve_all_recorded_spikes()

# Save results
suffix = end_time.strftime("_%H%M%S_%d%m%Y")
if args.filename:
    filename = args.filename
else:
    filename = "spinncer_experiment" + str(suffix)

# Check if the results folder exist
if not os.path.isdir(args.result_dir) and not os.path.exists(args.result_dir):
    os.mkdir(args.result_dir)

# Retrieve simulation parameters for provenance tracking and debugging purposes
sim_params = {
    "cell_params": cerebellum_circuit.retrieve_cell_params(),
    "argparser": vars(args),
    "git_hash": retrieve_git_commit(),
    "run_end_time": end_time.strftime("%H:%M:%S_%d/%m/%Y"),
    "wall_clock_run_time": str(total_time),
    "wall_clock_sim_run_time": str(sim_total_time),
}

# Save results to file in [by default] the `results/' directory
np.savez_compressed(os.path.join(args.result_dir, filename),
                    simulation_parameters=sim_params,
                    all_spikes=recorded_spikes,
                    stimulus_params=stimulus_params,
                    simtime=args.simtime,
                    **recorded_spikes)

# Appropriately end the simulation
sim.end()

# Report time taken
print("Results stored in  -- " + filename)

# Report time taken
print("Total time elapsed -- " + str(total_time))
