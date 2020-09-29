# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *
from spinncer.utilities.constants import *
from spinncer.utilities.neo_converter import convert_spikes
from spinncer.utilities import create_poisson_spikes, round_to_nearest_accum

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
# analysis functions
from spinncer.cerebellum_analysis import *
import pylab as plt
import sys
import os
import pandas as pd
import numpy as np
import traceback

# Record SCRIPT start time (wall clock)
start_time = plt.datetime.datetime.now()

# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=1,
          timescale=args.timescale)

if spinnaker_sim:
    sim.set_number_of_neurons_per_core(sim.IF_cond_exp, 5)

simtime = args.simtime

# 40 values for i_offset
dc_currents = np.linspace(0, 1., 61)[1:]
# use SpiNNaker representable values by casting to the closest accum
spinnaker_dc_currents = np.array(
    [round_to_nearest_accum(x, 0) for x in dc_currents]
)

# 10 values for sub-cycles
subcycles = np.arange(11)[1:]

n_neurons = dc_currents.size

typical_pop_dict = {k: None for k in CELL_PARAMS.keys() if k != "glomerulus"}

per_pop_r_mem = {}
pop_by_subcycle = {k: copy.deepcopy(typical_pop_dict) for k in subcycles}
recorded_spikes = copy.deepcopy(pop_by_subcycle)
recorded_voltage = copy.deepcopy(pop_by_subcycle)
cannonical_rbls = RMEM_RBLS if args.r_mem else VANILLA_RBLS

print("=" * 80)
print("\nThis is a test of the SpiNNaker solver when confronted with very rapidly changing sub-threshold voltage "
      "\nvalues to due high values on intrinsic current (DC, i_offset) injected into cells.")
print("Original DC current values", dc_currents, "nA")
print("SpiNNaker DC current values", spinnaker_dc_currents, "nA")
print("Diff in DC current values (Prescribed - S1615)", dc_currents-spinnaker_dc_currents)
print("Solver sub-cycle values", subcycles)
print("=" * 80)
for sc in subcycles:
    print("Creating populations using ", sc, "solver sub-cycles.")
    # Set up correct additional params
    additional_parameters = {
        "additional_parameters": {
            # No RB LS here as these connections don't receive any spikes
            "n_steps_per_timestep": sc
        }
    }
    # Set up all the pops for the current value of sub-cycles
    for pop in typical_pop_dict.keys():
        print("\tCreating pop", pop)
        curr_cell_params = copy.deepcopy(CELL_PARAMS[pop])
        # Do we have to consider the R_mem case?
        if args.r_mem:
            r_mem = curr_cell_params['tau_m'] / curr_cell_params['cm']
            per_pop_r_mem[pop] = r_mem
        else:
            per_pop_r_mem[pop] = 1

        # Set up correct cell parameters
        curr_cell_params['i_offset'] = spinnaker_dc_currents * per_pop_r_mem[pop]

        if spinnaker_sim:
            curr_pop = sim.Population(
                n_neurons,
                cellclass=sim.IF_cond_exp,
                cellparams=curr_cell_params,
                label="{} at {} sub-cycles".format(pop, sc),
                **additional_parameters)
        else:
            curr_pop = sim.Population(
                n_neurons,
                cellclass=sim.IF_cond_exp,
                cellparams=curr_cell_params,
                label="{} at {} sub-cycles".format(pop, sc))

        # Set up recordings here
        curr_pop.record(["spikes", "v"])

        # Add the current population to an object that might be a bit easier to parse
        pop_by_subcycle[sc][pop] = curr_pop

        # initialise V
        curr_pop.initialize(v=curr_cell_params['v_rest'])

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
for sc, pops_for_sc in pop_by_subcycle.items():
    for pop_label, pop_o in pops_for_sc.items():
        print("Retrieving spikes for ", pop_label, "...")
        recorded_spikes[sc][pop_label] = convert_spikes(pop_o.get_data(['spikes']))
        print("Retrieving v for ", pop_label, "...")
        recorded_voltage[sc][pop_label] = np.array(
            pop_o.get_data(['v']).segments[0].filter(name='v'))[0].T

sim.end()

# Check if the results folder exist
if not os.path.isdir(args.result_dir) and not os.path.exists(args.result_dir):
    os.mkdir(args.result_dir)

# save required csv
filename = "DC_testing_all_cells_{}_r_mem_{}".format(
    args.simulator,
    args.r_mem,
)
if args.suffix:
    filename += "_" + args.suffix

excel_filename = filename + ".xlsx"
np.savez_compressed(os.path.join(args.result_dir, filename),
                    # Recordings
                    spikes=recorded_spikes,
                    v=recorded_voltage,
                    # Experiment parameters
                    subcycles=subcycles,
                    simtime=simtime,
                    dc_currents=dc_currents,
                    spinnaker_dc_currents=spinnaker_dc_currents
                    )

print("R_mem for all pops")
pp(per_pop_r_mem)
print("Total time elapsed -- " + str(total_time))
