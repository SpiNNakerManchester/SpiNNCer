# argparser for easily running experiment from cli
from spinncer.spinncer_argparser import *
from spinncer.utilities.constants import *
from spinncer.utilities.neo_converter import convert_spikes
from spinncer.utilities import create_poisson_spikes, round_to_nearest_accum, floor_spike_time

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

SINGLE_SPIKE = {k: 1 for k in EXPECTED_MAX_SPIKES.keys()}

test_case_names = ["Single Spike", "Max Spike", "Full Spikes"]
cases = [0, 1]
test_spikes = [SINGLE_SPIKE, SINGLE_SPIKE, EXPECTED_MAX_SPIKES]
scale_weights = [SINGLE_SPIKE, EXPECTED_MAX_SPIKES, SINGLE_SPIKE]

per_pop_incoming_projections = {}
for k in CELL_PARAMS.keys():
    per_pop_incoming_projections[k] = []

for proj, proj_params in CONNECTIVITY_MAP.items():
    post_pop = proj_params['post']
    per_pop_incoming_projections[post_pop].append(proj)

# Set up the simulation
sim.setup(timestep=args.timestep, min_delay=args.timestep, max_delay=1,
          timescale=args.timescale,
          spike_precision=args.nest_grid  # NEST Spike precision
          )

if spinnaker_sim:
    sim.set_number_of_neurons_per_core(sim.IF_cond_exp, 5)

simtime = args.simtime

# 4 values for sub-cycles. If we really need more we'll find out and modify the value below.
subcycles = np.arange(11)[1:]

n_neurons = 1
n_test_cells = 1

spike_time = 10 - args.timestep

spike_source_per_no_spikes = {}

typical_pop_dict = {k: None for k in CONNECTIVITY_MAP.keys() if k != "glomerulus"}
cases_dict = {k: copy.deepcopy(typical_pop_dict) for k in cases}

per_pop_r_mem = {}
pop_by_subcycle = {k: copy.deepcopy(cases_dict) for k in subcycles}

input_spikes_by_case = {k: {} for k in cases}

recorded_spikes = copy.deepcopy(pop_by_subcycle)
recorded_voltage = copy.deepcopy(pop_by_subcycle)
recorded_gsyn_exc = copy.deepcopy(pop_by_subcycle)
recorded_gsyn_inh = copy.deepcopy(pop_by_subcycle)

cannonical_rbls = RMEM_RBLS if args.r_mem else VANILLA_RBLS

pop_order = get_plot_order(EXPECTED_MAX_SPIKES.keys())

for case, test_name, spikes_for_test, do_weight_scaling in \
        zip(cases, test_case_names, test_spikes, scale_weights):
    print("Current case:", test_name)
    for sc in subcycles:
        print("Creating populations using ", sc, "solver sub-cycles.")
        # Set up all the pops for the current value of sub-cycles
        # seed the elephant spike generator here so that all sub-cycles generate the same spikes
        np.random.seed(3141592653)
        # Consistently loop over populations in the same order
        for proj in pop_order:
            # Set up correct additional params
            pop = CONNECTIVITY_MAP[proj]['post']
            rb_ls = cannonical_rbls[pop]
            additional_parameters = {
                "additional_parameters": {
                    "rb_left_shifts": rb_ls,
                    "n_steps_per_timestep": sc
                }
            }

            print("\tCreating pop", proj)
            curr_cell_params = copy.deepcopy(CELL_PARAMS[pop])
            # Do we have to consider the R_mem case?
            if args.r_mem and spinnaker_sim:
                r_mem = curr_cell_params['tau_m'] / curr_cell_params['cm']
                per_pop_r_mem[proj] = r_mem
            else:
                per_pop_r_mem[proj] = 1

            # Set up correct cell parameters
            if 'i_offset' in curr_cell_params.keys():
                curr_cell_params['i_offset'] = round_to_nearest_accum(
                    curr_cell_params['i_offset'] * per_pop_r_mem[proj], 0)

            if spinnaker_sim:
                curr_pop = sim.Population(
                    n_test_cells,  # Testing a single cell
                    cellclass=sim.IF_cond_exp,
                    cellparams=curr_cell_params,
                    label="{}-{} at {} sub-cycles".format(test_name, proj, sc),
                    **additional_parameters)
            else:
                curr_pop = sim.Population(
                    n_test_cells,  # Testing a single cell
                    cellclass=sim.IF_cond_exp,
                    cellparams=curr_cell_params,
                    label="{}-{} at {} sub-cycles".format(test_name, proj, sc))

            # Set up recordings here
            curr_pop.record(["spikes", "v", "gsyn_inh", "gsyn_exc"])

            # Add the current population to an object that might be a bit easier to parse
            pop_by_subcycle[sc][case][proj] = curr_pop

            # initialise V
            curr_pop.initialize(v=curr_cell_params['v_rest'])

            # Retrieve population level firing rates
            no_spikes = spikes_for_test[proj]

            curr_weight = CONNECTIVITY_MAP[proj]['weight']
            is_inh = int(curr_weight < 0)
            curr_weight = np.abs(curr_weight) * per_pop_r_mem[proj]
            # Extra scaling for weights for case 1 (Single spikes, max weight)
            curr_weight *= do_weight_scaling[proj]

            if no_spikes not in spike_source_per_no_spikes.keys():
                # add the exc and inh input populations
                ssa = sim.Population(
                    n_neurons,
                    cellclass=sim.SpikeSourceArray,
                    cellparams={
                        'spike_times': [spike_time] * no_spikes
                    },
                    label="{}-SSA for {} at {} sub-cycles".format(test_name, proj, sc))
                spike_source_per_no_spikes[no_spikes] = ssa
            else:
                ssa = spike_source_per_no_spikes[no_spikes]

            input_spikes_by_case[case][proj] = [spike_time] * no_spikes

            # add the exc and inh projections
            sim.Projection(ssa, curr_pop,
                           sim.AllToAllConnector(),
                           synapse_type=sim.StaticSynapse(weight=curr_weight),
                           receptor_type='inhibitory' if is_inh else 'excitatory',
                           label="{}-SSA projection to {} at {} sub-cycles".format(test_name, proj, sc)
                           )
    print("-" * 80)

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
    for case, pops_for_case in pops_for_sc.items():
        for pop_label, pop_o in pops_for_case.items():
            print("Retrieving spikes for ", pop_label, "...")
            recorded_spikes[sc][case][pop_label] = convert_spikes(pop_o.get_data(['spikes']))

            print("Retrieving v for ", pop_label, "...")
            recorded_voltage[sc][case][pop_label] = np.array(
                pop_o.get_data(['v']).segments[0].filter(name='v'))[0].T

            print("Retrieving gsyn exc for ", pop_label, "...")
            recorded_gsyn_exc[sc][case][pop_label] = np.array(
                pop_o.get_data(['gsyn_exc']).segments[0].filter(name='gsyn_exc'))[0].T / per_pop_r_mem[pop_label]

            print("Retrieving gsyn inh for ", pop_label, "...")
            recorded_gsyn_inh[sc][case][pop_label] = np.array(
                pop_o.get_data(['gsyn_inh']).segments[0].filter(name='gsyn_inh'))[0].T / per_pop_r_mem[pop_label]

sim.end()

# Check if the results folder exist
if not os.path.isdir(args.result_dir) and not os.path.exists(args.result_dir):
    os.mkdir(args.result_dir)

# save required csv
filename = "SINGLE_SPIKE_testing_all_cells_{}_r_mem_{}".format(
    args.simulator,
    args.r_mem,
)
if args.suffix:
    filename += "_" + args.suffix

excel_filename = filename + ".xlsx"
np.savez_compressed(os.path.join(args.result_dir, filename),
                    # Recordings
                    # Spikes
                    spikes=recorded_spikes,

                    # Others
                    v=recorded_voltage,
                    gsyn_exc=recorded_gsyn_exc,
                    gsyn_inh=recorded_gsyn_inh,

                    # Experiment parameters
                    subcycles=subcycles,
                    simtime=simtime,
                    input_spikes=input_spikes_by_case,
                    cases=cases,
                    case_names=test_case_names
                    )

print("R_mem for all pops")
pp(per_pop_r_mem)
print("Total time elapsed -- " + str(total_time))
