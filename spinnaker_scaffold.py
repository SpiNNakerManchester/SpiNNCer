try:
    import pyNN.spinnaker as sim
except:
    import spynnaker8 as sim

from copy import deepcopy

import h5py
import matplotlib.pyplot as plt
from scaffold_params import *


def connect_neuron(conn_mat, pre, post, syn_param):
    WEIGHT = syn_param["weight"]
    DELAY = syn_param["delay"]

    switcher = {'granule': 1117, 'golgi': 2, 'purkinje': 1119,
                'basket': 1121, 'stellate': 1123,
                'glomerulus': 1150, 'dcn': 1125}
    pre_idx = []
    post_idx = []

    for x in conn_mat[:, 0] + 1:
        fl = np.where(pre.all_cells - switcher[syn_param["pre"]] == x)

        if np.shape(fl[0])[0] != 0:
            pre_idx.append([fl[0][0]])

    for x in conn_mat[:, 1] + 1:
        fl = np.where(post.all_cells - switcher[syn_param["post"]] == x)

        post_idx.append(fl[0][0])

    conn_list = []
    for i in range(0, len(pre_idx) - 1):
        p = (pre_idx[i][0], post_idx[i], WEIGHT, DELAY)
        conn_list.append(p)

    conn = sim.FromListConnector(conn_list)
    ss = sim.StaticSynapse(weight=WEIGHT, delay=DELAY)

    sim.Projection(pre, post, conn, ss)


DELAY = 0.2
sim.setup(timestep=0.1, min_delay=0.1, max_delay=1.5)

filename = 'scaffold_detailed__158.0x158.0_v3.hdf5'
# filename = 'scaffold_full_dcn_400.0x400.0_v3.hdf5'
f = h5py.File(filename, 'r+')
positions = np.array(f['positions'])
num_tot = 0

sorted_nrn_types = sorted(list(cell_type_ID.values()))

# Create a 'inverse map' of cell_type_IDaa_pc
id_2_cell_type = {val: key for key, val in cell_type_ID.items()}

# Create a dictionary with all cell names (keys)
# and lists that will contain nest models (values)
neuron_models = {key: [] for key in cell_type_ID.keys()}
no_neurons = {key: [] for key in cell_type_ID.keys()}
per_model_cell_params = {key: [] for key in cell_type_ID.keys()}

for cell_id in sorted_nrn_types:
    cell_name = id_2_cell_type[cell_id]
    if cell_name != 'glomerulus':

        if cell_name == 'golgi':
            cell_params = {'tau_refrac': 2.0,  # ms
                           'cm': 0.076,  # nF
                           'v_thresh': -55.0,  # mV
                           'v_reset': -75.0,  # mV
                           'tau_m': 21.1,
                           # 'e_rev_leak': -65.0,  # mV
                           'i_offset': 36.75,  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
                           'tau_syn_E': 0.5,
                           'tau_syn_I': 10.0}
        elif cell_name == 'granule':
            cell_params = {'tau_refrac': 1.5,  # ms
                           'cm': 0.003,  # nF
                           'v_thresh': -42.0,  # mV
                           'v_reset': -84.0,  # mV
                           'tau_m': 20,
                           # 'v_rest': -74.0,  # mV
                           'i_offset': 0,  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
                           'tau_syn_E': 0.5,
                           'tau_syn_I': 10.0}

        elif cell_name == 'purkinje':
            cell_params = {'tau_refrac': 0.8,  # ms
                           'cm': 0.062,  # nF
                           'v_thresh': -47.0,  # mV
                           'v_reset': -72.0,  # mV
                           'tau_m': 88.6,
                           # 'e_rev_leak': -62.0,  # mV
                           'i_offset': 750,  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
                           'tau_syn_E': 0.5,
                           'tau_syn_I': 10.0}
        elif cell_name == 'stellate' or cell_name == 'basket':
            cell_params = {'tau_refrac': 1.59,  # ms
                           'cm': 0.0146,  # nF
                           'v_thresh': -53.0,  # mV
                           'v_reset': -78.0,  # mV
                           'tau_m': 14.6,
                           # 'e_rev_leak': -68.0,  # mV
                           'i_offset': 15.6,  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
                           'tau_syn_E': 0.5,
                           'tau_syn_I': 10.0}
        elif cell_name == 'dcn':
            cell_params = {'tau_refrac': 3.7,  # ms
                           'cm': 0.089,  # nF
                           'v_thresh': -48.0,  # mV
                           'v_reset': -69.0,  # mV
                           'tau_m': 57.1,
                           # 'e_rev_leak': -59.0,  # mV
                           'i_offset': 45.75,  # pA # tonic ~9-10 Hz  ;previous = 36.0 pA
                           'tau_syn_E': 0.5,
                           'tau_syn_I': 10.0}

        model = sim.IF_cond_exp(**cell_params)

    else:
        # model = SpikeSourcePoisson()
        model = sim.IF_cond_exp(**cell_params)

    cell_pos = positions[positions[:, 1] == cell_id, :]
    neuron_models[cell_name] = sim.Population(cell_pos.shape[0], model)
    # Recording useful values
    no_neurons[cell_name] = cell_pos.shape[0]
    per_model_cell_params[cell_name] = deepcopy(cell_params)
    if cell_name == "granule":
        pop_gr = neuron_models[cell_name]
    num_tot += cell_pos.shape[0]

# noise = NoisyCurrentSource(mean=0, stdev=50)
# print(neuron_models.get('dcn'))
# print(neuron_models['dcn'].all_cells)
# print(neuron_models['dcn'][1])


# for name in neuron_models.keys():
#    a = neuron_models.get(name)
#    a.inject(noise)
#    noise.inject_into(cell)


conn_aa_goc = np.array(f['connections/aa_goc'])
syn_p = {"model": "static_synapse", "weight": 20, "delay": 0.5, "pre": "granule", "post": "golgi"}
# connect_neuron(conn_aa_goc, neuron_models.get('granule'), neuron_models.get('golgi'), syn_p)

conn_aa_pc = np.array(f['connections/aa_pc'])
syn_p = {"model": "static_synapse", "weight": 50.0, "delay": 0.5, "pre": "granule", "post": "purkinje"}
# connect_neuron(conn_aa_pc, neuron_models.get('granule'), neuron_models.get('purkinje'), syn_p)

conn_bc_pc = np.array(f['connections/bc_pc'])
syn_p = {"model": "static_synapse", "weight": -2.50, "delay": 0.5, "pre": "basket", "post": "purkinje"}
# connect_neuron(conn_bc_pc, neuron_models.get('basket'), neuron_models.get('purkinje'), syn_p)

conn_gj_bc = np.array(f['connections/gj_bc'])
syn_p = {"model": "static_synapse", "weight": -3, "delay": 0.5, "pre": "basket", "post": "basket"}
# connect_neuron(conn_gj_bc, neuron_models.get('basket'), neuron_models.get('basket'), syn_p)

conn_gj_sc = np.array(f['connections/gj_sc'])
syn_p = {"model": "static_synapse", "weight": -2.5, "delay": 0.5, "pre": "stellate", "post": "stellate"}
# connect_neuron(conn_gj_sc, neuron_models.get('stellate'), neuron_models.get('stellate'), syn_p)

conn_glom_goc = np.array(f['connections/glom_goc'])
syn_p = {"model": "static_synapse", "weight": 2.0, "delay": 0.5, "pre": "glomerulus", "post": "golgi"}
# connect_neuron(conn_glom_goc, neuron_models.get('glomerulus'), neuron_models.get('golgi'), syn_p)
# syn_p = {"model": "static_synapse", "weight": 2.0, "delay": 0.5, "pre": "0", "post": "golgi"}
# connect_neuron(conn_glom_goc, stimulus, neuron_models.get('golgi'), syn_p)

conn_glom_grc = np.array(f['connections/glom_grc'])
syn_p = {"model": "static_synapse", "weight": 9.0, "delay": 0.5, "pre": "glomerulus", "post": "granule"}
# connect_neuron(conn_glom_grc, neuron_models.get('glomerulus'), neuron_models.get('granule'), syn_p)

conn_goc_grc = np.array(f['connections/goc_grc'])
syn_p = {"model": "static_synapse", "weight": -5, "delay": 0.5, "pre": "golgi", "post": "granule"}
# connect_neuron(conn_goc_grc, neuron_models.get('golgi'), neuron_models.get('granule'), syn_p)

conn_pc_dcn = np.array(f['connections/pc_dcn'])
syn_p = {"model": "static_synapse", "weight": -0.2, "delay": 0.5, "pre": "purkinje", "post": "dcn"}
# connect_neuron(conn_pc_dcn, neuron_models.get('purkinje'), neuron_models.get('dcn'), syn_p)

conn_pf_bc = np.array(f['connections/pf_bc'])
syn_p = {"model": "static_synapse", "weight": 0.4, "delay": 1, "pre": "granule", "post": "basket"}
# connect_neuron(conn_pf_bc, neuron_models.get('granule'), neuron_models.get('basket'), syn_p)

conn_pf_goc = np.array(f['connections/pf_goc'])
syn_p = {"model": "static_synapse", "weight": 0.4, "delay": 1, "pre": "granule", "post": "golgi"}
# connect_neuron(conn_pf_goc, neuron_models.get('granule'), neuron_models.get('golgi'), syn_p)

conn_pf_pc = np.array(f['connections/pf_pc'])
syn_p = {"model": "static_synapse", "weight": 0.05, "delay": 1, "pre": "granule", "post": "purkinje"}
# connect_neuron(conn_pf_pc, neuron_models.get('granule'), neuron_models.get('purkinje'), syn_p)

conn_pf_sc = np.array(f['connections/pf_sc'])
syn_p = {"model": "static_synapse", "weight": 0.4, "delay": 1, "pre": "granule", "post": "stellate"}
# connect_neuron(conn_pf_sc, neuron_models.get('granule'), neuron_models.get('stellate'), syn_p)

conn_sc_pc = np.array(f['connections/sc_pc'])
syn_p = {"model": "static_synapse", "weight": -2, "delay": 0.5, "pre": "stellate", "post": "purkinje"}
# connect_neuron(conn_sc_pc, neuron_models.get('stellate'), neuron_models.get('purkinje'), syn_p)

conn_glom_dcn = np.array(f['connections/glom_dcn'])
syn_p = {"model": "static_synapse", "weight": 2.0, "delay": 0.5, "pre": "glomerulus", "post": "dcn"}
# connect_neuron(conn_glom_dcn, neuron_models.get('glomerulus'), neuron_models.get('dcn'), syn_p)

conn_gj_goc = np.array(f['connections/gj_goc'])
syn_p = {"model": "static_synapse", "weight": -8.0, "delay": 0.5, "pre": "golgi", "post": "golgi"}
# connect_neuron(conn_gj_goc, neuron_models.get('golgi'), neuron_models.get('golgi'), syn_p)

###########  Stim and simulation

TOT_DURATION = 500.  # mseconds
STIM_START = 300.  # beginning of stimulation
STIM_END = 350.  # end of stimulation
STIM_FREQ = 120.  # Frequency in Hz
RADIUS = 50.5  # Microm

# setup(timestep= 1.0)
spike_nums = np.int(np.round((STIM_FREQ * (STIM_END - STIM_START)) / TOT_DURATION))
# spike_nums = len(neuron_models.get('golgi'))
stim_array = np.round(np.linspace(STIM_START, STIM_END, spike_nums))
stim_array_int = []
for i in stim_array:
    stim_array_int.append(int(i))
origin = np.array([200., 200., 75.])
phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0 * np.pi:100j]
x = RADIUS * np.sin(phi) * np.cos(theta) + origin[0]
y = RADIUS * np.sin(phi) * np.sin(theta) + origin[1]
z = RADIUS * np.cos(phi) + origin[2]

gloms_pos = positions[positions[:, 1] == cell_type_ID['glomerulus'], :]
# find center of 'glomerular sphere'
x_c, y_c, z_c = np.median(gloms_pos[:, 2]), np.median(gloms_pos[:, 3]), np.median(gloms_pos[:, 4])

# Find glomeruli falling into the selected volume
target_gloms_idx = np.sum((gloms_pos[:, 2::] - np.array([x_c, y_c, z_c])) ** 2, axis=1).__lt__(RADIUS ** 2)
tg_idx = []
# for i in range(neuron_models.get('granule').size): #range(len(target_gloms_idx)):
#    if target_gloms_idx[i] == True:
#        tg_idx.append(i)

# target_gloms = gloms_pos[tg_idx, 0] + 1
# print(neuron_models['glomerulus'].all_cells)
# id_stim = [glom for glom in neuron_models['glomerulus'].all_cells if glom in target_gloms]

# count = 0
# id_stim = []
# for glom in neuron_models['glomerulus'].all_cells:
#    print(glom in target_gloms)
#    if(glom in target_gloms):
#        id_stim.append(count)
#    for i in target_gloms:
#        if(i == glom):
#            id_stim.append(count)
#    count += 1


# n = len(target_gloms)

print(neuron_models.get('granule').size)
# spikeTimes= [2, 4, 6, 8, 10]
stimulus = sim.Population(pop_gr.size, sim.SpikeSourceArray(spike_times=stim_array_int))
# neuron_models.get('granule').size
conn_list = []

# for i in range(n):
#    p = (0, tg_idx[i], 1., 1)
#    conn_list.append(p)

# conn = ListConnector(conn_list)
ss = sim.StaticSynapse(weight=1)
# Projection(stimulus, neuron_models.get('glomerulus'), conn, ss)
# Projection(stimulus, neuron_models.get('granule'), conn, ss)
# Projection(stimulus, neuron_models.get('glomerulus')[tg_idx], AllToAllConnector(allow_self_connections=False), ss)

# stimulus = Population(1, SpikeSourceArray(spike_times= stim_array_int))
sim.Projection(stimulus, pop_gr, sim.OneToOneConnector(), ss)

# for i in neuron_models.values():
#    i.record(['v'])
# pop_gr.record(['v'])

# RECORD SPIKES FOR ALL POPULATIONS
for cell_id in sorted_nrn_types:
    cell_name = id_2_cell_type[cell_id]
    print("Recording", cell_name, "...")
    neuron_models[cell_name].record(['spikes'])


sim.run(TOT_DURATION)

# g_data = neuron_models.get('granule')
# data1 = pop_gr.get_data(variables=["v"])

# Figure(Panel(g_data.get_data('spikes')), title="Prova")

# print(data1)

recorded_spikes = {}
for cell_id in sorted_nrn_types:
    cell_name = id_2_cell_type[cell_id]
    print("Retrieving recording for", cell_name, "...")
    recorded_spikes[cell_name] = neuron_models[cell_name].spinnaker_get_data('spikes')

np.savez_compressed("results_for_scaffold_experiment",
                    spikes=recorded_spikes,
                    network_filename=filename,
                    simtime=TOT_DURATION,
                    no_neurons=no_neurons,
                    per_model_cell_params=per_model_cell_params)
