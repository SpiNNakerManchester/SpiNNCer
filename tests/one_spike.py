"""
This script tests whether spike counting additional provenance is correct
"""
import pyNN.spiNNaker as sim
sim.setup(timestep=1, min_delay=1, max_delay=1)

# Compute 1 spike
n_neurons = 100
spike_times = [[] for _ in range(n_neurons)]
spike_times[42].append(314)

# input pop
inp_pop = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': spike_times},
    label="One spike pop")

# lif pop
lif_pop = sim.Population(
    100,
    cellclass=sim.IF_curr_exp,
    label="LIF pop 1 spike")

# connection
sim.Projection(inp_pop, lif_pop, sim.OneToOneConnector(),
               synapse_type=sim.StaticSynapse(weight=1.0, delay=1.0))

sim.run(1000)
sim.end()
