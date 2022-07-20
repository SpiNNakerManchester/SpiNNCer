"""
This script tests whether delays are implemented correctly
"""
import pyNN.spynnaker as sim
from pyNN.utility.plotting import Figure, Panel
import matplotlib.pyplot as plt
import numpy as np

sim.setup(timestep=0.1, min_delay=0.1, max_delay=80)

# Compute 100 spike
n_neurons = 5
runtime = 50  # ms
spike_times = [[0] for _ in range(n_neurons)]
DELAY = 5
WEIGHT = 5
# Compute 1:1 connectivity
conn = []
for i in range(n_neurons):
    conn.append((i, i))

# input pop
inp_pop = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': spike_times},
    label="100 spike pop")

# lif pop
lif_pop = sim.Population(
    n_neurons,
    cellclass=sim.IF_curr_exp,
    label="POP RECEIVING DELAYED STIMULUS")

lif_pop.record(['v', 'spikes'])
# connection
proj = sim.Projection(inp_pop, lif_pop, sim.FromListConnector(conn),
                      synapse_type=sim.StaticSynapse(
                          weight=WEIGHT, delay=DELAY))

sim.run(runtime)

# get data (could be done as one, but can be done bit by bit as well)
v = lif_pop.get_data('v')
spikes = lif_pop.get_data('spikes')
retrieved_conn = np.array(proj.get(('weight', 'delay'),
                                   format="list")._get_data_items())
Figure(
    # raster plot of the presynaptic neuron spike times
    Panel(spikes.segments[0].spiketrains,
          yticks=True, markersize=0.5, xlim=(0, runtime)),
    # membrane potential of the postsynaptic neuron
    Panel(v.segments[0].filter(name='v')[0],
          ylabel="Membrane potential (mV)",
          data_labels=[lif_pop.label], yticks=True, xlim=(0, runtime),
          xlabel="Time (ms)", xticks=True),
    title="Simple synfire chain example",
    annotations="Simulated with {}".format(sim.name())
)
plt.show()

sim.end()

assert all(retrieved_conn['delay'] == DELAY), retrieved_conn['delay']
assert all(retrieved_conn['weight'] == WEIGHT), retrieved_conn['weight']
