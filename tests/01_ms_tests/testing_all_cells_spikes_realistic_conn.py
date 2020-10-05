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

# Information on firing rates used for each pop during tests
low_case_rates = {
    'glomerulus': 0.92,
    'granule': 2.05,
    'golgi': 18.63,
    'stellate': 31.68,
    'basket': 27.93,
    'purkinje': 47.68,
    'dcn': 17.74
}

high_case_rates = {
    'glomerulus': 144.38,
    'granule': 89.59,
    'golgi': 135.13,
    'stellate': 220.75,
    'basket': 193.04,
    'purkinje': 381.82,
    'dcn': 0
}

test_rates = [low_case_rates, high_case_rates]
test_case_names = ["LOW", "HIGH"]
cases = [0, 1]

# need information on peak fan-in from each connection.
max_fan_ins = {'aa_goc': 400,
 'aa_pc': 278,
 'bc_pc': 20,
 'gj_bc': 12,
 'gj_goc': 50,
 'gj_sc': 11,
 'glom_dcn': 147,
 'glom_goc': 110,
 'glom_grc': 4,
 'goc_grc': 4,
 'pc_dcn': 30,
 'pf_bc': 1343,
 'pf_goc': 1600,
 'pf_pc': 29196,
 'pf_sc': 1379,
 'sc_pc': 20}


per_pop_incoming_projections = {}
for k in low_case_rates.keys():
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

n_test_cells = 1

typical_pop_dict = {k: None for k in CELL_PARAMS.keys() if k != "glomerulus"}
cases_dict = {k: copy.deepcopy(typical_pop_dict) for k in cases}

spike_source_per_pop = {k: {} for k in cases}

per_pop_r_mem = {}
pop_by_subcycle = {k: copy.deepcopy(cases_dict) for k in subcycles}

input_spikes_by_case = {k: {} for k in cases}

recorded_spikes = copy.deepcopy(pop_by_subcycle)
recorded_voltage = copy.deepcopy(pop_by_subcycle)
recorded_gsyn_exc = copy.deepcopy(pop_by_subcycle)
recorded_gsyn_inh = copy.deepcopy(pop_by_subcycle)

cannonical_rbls = RMEM_RBLS if args.r_mem else VANILLA_RBLS

print("=" * 80)
print(
    "\nThis is a test of the SpiNNaker solver when confronted with activity evoking realistic peak synaptic conductance.")
print("LOW population firing rates", low_case_rates, "Hz")
print("HIGH population firing rates", high_case_rates, "Hz")
print("Solver sub-cycle values", subcycles)
print("=" * 80)

pop_order = get_plot_order(low_case_rates.keys())

# Set up all the pops for the current value of sub-cycles
# seed the elephant spike generator here so that all sub-cycles generate the same spikes
np.random.seed(3141592653)

for case, test_name, rates_for_test in zip(cases, test_case_names, test_rates):
    print("Current case:", test_name)
    for sc in subcycles:
        print("Creating populations using ", sc, "solver sub-cycles.")
        # Consistently loop over populations in the same order
        for pop in pop_order:
            # Set up correct additional params
            if pop in ['glomerulus']:
                continue
            rb_ls = cannonical_rbls[pop]
            additional_parameters = {
                "additional_parameters": {
                    "rb_left_shifts": rb_ls,
                    "n_steps_per_timestep": sc
                }
            }

            print("\tCreating pop", pop)
            curr_cell_params = copy.deepcopy(CELL_PARAMS[pop])
            # Do we have to consider the R_mem case?
            if args.r_mem and spinnaker_sim:
                r_mem = curr_cell_params['tau_m'] / curr_cell_params['cm']
                per_pop_r_mem[pop] = r_mem
            else:
                per_pop_r_mem[pop] = 1

            # Set up correct cell parameters
            if 'i_offset' in curr_cell_params.keys():
                curr_cell_params['i_offset'] = round_to_nearest_accum(
                    curr_cell_params['i_offset'] * per_pop_r_mem[pop], 0)

            if spinnaker_sim:
                curr_pop = sim.Population(
                    n_test_cells,  # Testing a single cell
                    cellclass=sim.IF_cond_exp,
                    cellparams=curr_cell_params,
                    label="{}-{} at {} sub-cycles".format(test_name, pop, sc),
                    **additional_parameters)
            else:
                curr_pop = sim.Population(
                    n_test_cells,  # Testing a single cell
                    cellclass=sim.IF_cond_exp,
                    cellparams=curr_cell_params,
                    label="{}-{} at {} sub-cycles".format(test_name, pop, sc))

            # Set up recordings here
            curr_pop.record(["spikes", "v", "gsyn_inh", "gsyn_exc"])

            # Add the current population to an object that might be a bit easier to parse
            pop_by_subcycle[sc][case][pop] = curr_pop

            # initialise V
            curr_pop.initialize(v=curr_cell_params['v_rest'])

            # retrieve the proportions of each weight
            inc_proj = per_pop_incoming_projections[pop]
            inc_proj.sort()

            # Create a spike source array for each pre-synaptic pop
            for proj in inc_proj:
                curr_weight = CONNECTIVITY_MAP[proj]['weight']
                is_inh = int(curr_weight < 0)
                curr_weight = np.abs(curr_weight) * per_pop_r_mem[pop]
                pre_pop = CONNECTIVITY_MAP[proj]['pre']
                # Retrieve population level firing rates
                rate_for_pre_pop = rates_for_test[pre_pop]

                # Go from population-level rate to individual neuron rates
                curr_fan_in = int(max_fan_ins[proj])

                if curr_fan_in > 100:
                    print("Current fan in is", curr_fan_in, "which is rather high")
                    curr_fan_in = int(curr_fan_in / 10.)
                    rate_for_pre_pop *= 10

                if pre_pop not in spike_source_per_pop[case].keys():
                    # Produce spikes
                    curr_spikes = create_poisson_spikes(curr_fan_in,
                                                        [[rate_for_pre_pop], ] * curr_fan_in,
                                                        [[0], ] * curr_fan_in,
                                                        [[simtime], ] * curr_fan_in)

                    # Round spike times to time step boundary
                    for id, exc_s in enumerate(curr_spikes):
                        rounded_spike_times = floor_spike_time(exc_s, t_stop=simtime * pq.ms)
                        # DEALING WITH nest.lib.hl_api_exceptions.NESTErrors.BadProperty:
                        # ("BadProperty in SetStatus_id: Setting status of a
                        # 'spike_generator' with GID 855: spike time cannot be
                        # set to 0.", 'BadProperty',
                        # <SLILiteral: SetStatus_id>, ": Setting status of a
                        # 'spike_generator' with GID 855: spike time cannot be set to 0.")
                        # Which means IT CAN'T BE 0.1, NOT 0
                        rounded_spike_times[rounded_spike_times < 0.2] = 0.2
                        curr_spikes[id] = rounded_spike_times

                    # save the spikes but don't record them from the Spike Sources
                    input_spikes_by_case[case][proj] = curr_spikes

                    # add the exc and inh input populations
                    ssa = sim.Population(
                        curr_fan_in,
                        cellclass=sim.SpikeSourceArray,
                        cellparams={
                            'spike_times': curr_spikes
                        },
                        label="{}-SSA for {} at {} sub-cycles".format(test_name, proj, sc))

                    spike_source_per_pop[case][pre_pop] = ssa
                else:
                    ssa = spike_source_per_pop[case][pre_pop]
                # add the exc and inh projections
                sim.Projection(ssa, curr_pop,
                               sim.AllToAllConnector(),
                               synapse_type=sim.StaticSynapse(weight=curr_weight),
                               receptor_type='inhibitory' if is_inh else 'excitatory',
                               label="{}- SSA projection to {} at {} sub-cycles".format(test_name, proj, sc)
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
filename = "SPIKES_REAL_FAN_INS_testing_all_cells_{}_r_mem_{}".format(
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
