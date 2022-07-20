"""
This script tests whether spike counting additional provenance is correct
"""
import pyNN.spynnaker as sim
sim.setup(timestep=1, min_delay=1, max_delay=1)


# input pop
inp_pop = sim.Population(
    100,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': []},
    label="No spikes pop")

# lif pop
lif_pop = sim.Population(
    100,
    cellclass=sim.IF_curr_exp,
    label="LIF pop 0 spikes")

# connection
sim.Projection(inp_pop, lif_pop, sim.OneToOneConnector(),
               synapse_type=sim.StaticSynapse(weight=1.0, delay=1.0))

sim.run(1000)
sim.end()
