"""
This script tests whether spike counting additional provenance is correct
"""
import pyNN.spiNNaker as sim
sim.setup(timestep=0.1, min_delay=0.1, max_delay=1)

# Generate 6000 spikes (1 for each input neuron)
n_neurons = 6000
spike_times = [[314] for _ in range(n_neurons)]
# Generate connectivity
conn = []
for i in range(n_neurons):
    conn.append((i, 0))

# input pop
inp_pop = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': spike_times},
    label="6000 spike pop")

# lif pop
lif_pop = sim.Population(
    1,
    cellclass=sim.IF_curr_exp,
    label="LIF pop")

# connection
sim.Projection(inp_pop, lif_pop, sim.FromListConnector(conn),
               synapse_type=sim.StaticSynapse(weight=1.0, delay=1.0))

sim.run(1000)
sim.end()
