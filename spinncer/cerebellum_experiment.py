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

# Record simulation start time (wall clock)
start_time = plt.datetime.datetime.now()
connectivity_filename = 'datasets/scaffold_detailed__158.0x158.0_v3.hdf5'

# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=10)

cerebellum_circuit = Cerebellum(sim, connectivity_filename)

# Test various exposed methods
populations = cerebellum_circuit.get_all_populations()
assert (len(populations) == 7)
input_populations = cerebellum_circuit.get_circuit_inputs()
assert (len(input_populations) == 1)
output_populations = cerebellum_circuit.get_circuit_outputs()
assert (len(output_populations) == 1)

# Add stimulus to the network

# Set up recordings
cerebellum_circuit.record_all_spikes()

# Run the simulation
sim.run(args.simtime)

# Retrieve recordings
recorded_spikes = cerebellum_circuit.retrieve_all_recorded_spikes()

# Compute time taken to reach this point
end_time = plt.datetime.datetime.now()
total_time = end_time - start_time

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
    "wall_clock_run_time": str(total_time)
}

# Save results to file in [by default] the `results/' directory
np.savez_compressed(os.path.join(args.result_dir, filename),
                    simulation_parameters=sim_params,
                    all_spikes=recorded_spikes,
                    simtime=args.simtime,
                    **recorded_spikes)

# Appropriately end the simulation
sim.end()

# Report time taken
print("Results stored in  -- " + filename)

# Report time taken
print("Total time elapsed -- " + str(total_time))
