# Copyright (c) 2017-2019 The University of Manchester
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Synfirechain-like example
"""
import pyNN.spiNNaker as sim
from pyNN.utility.plotting import Figure, Panel
import matplotlib.pyplot as plt

runtime = 5000
sim.setup(timestep=1.0, min_delay=1.0, max_delay=144.0)
nNeurons = 200  # number of neurons in each population
neuron_model = sim.IF_cond_alpha
sim.set_number_of_neurons_per_core(sim.IF_cond_alpha, nNeurons / 2)

cell_params_lif = {'cm': 0.25,
                   'i_offset': 0.0,
                   'tau_m': 20.0,
                   'tau_refrac': 2.0,
                   'tau_syn_E': 5.0,
                   'tau_syn_I': 5.0,
                   'v_reset': -70.0,
                   'v_rest': -65.0,
                   'v_thresh': -50.0,
                   'e_rev_E': 0.,
                   'e_rev_I': -80.
                   }

weight_to_spike = 0.2  # c.f. 0.035 for synfire_if_cond_exp.py
delay = 10

loopConnections = list()
for i in range(0, nNeurons):
    singleConnection = ((i, (i + 1) % nNeurons, weight_to_spike, delay))
    loopConnections.append(singleConnection)

injectionConnection = [(0, 0)]
spikeArray = {'spike_times': [[0]]}
main_pop = sim.Population(
    nNeurons, neuron_model(**cell_params_lif), label='pop_1')
input_pop = sim.Population(
    1, sim.SpikeSourceArray(**spikeArray), label='inputSpikes_1')

sim.Projection(
    main_pop, main_pop, sim.FromListConnector(loopConnections),
    sim.StaticSynapse(weight=weight_to_spike, delay=delay))
sim.Projection(
    input_pop, main_pop, sim.FromListConnector(injectionConnection),
    sim.StaticSynapse(weight=weight_to_spike, delay=1))

main_pop.record(['v', 'gsyn_exc', 'gsyn_inh', 'spikes'])

sim.run(runtime)

# get data (could be done as one, but can be done bit by bit as well)
v = main_pop.get_data('v')
gsyn_exc = main_pop.get_data('gsyn_exc')
gsyn_inh = main_pop.get_data('gsyn_inh')
spikes = main_pop.get_data('spikes')

figure_filename = "results.png"
Figure(
    # raster plot of the presynaptic neuron spike times
    Panel(spikes.segments[0].spiketrains,
          yticks=True, markersize=0.2, xlim=(0, runtime)),
    # membrane potential of the postsynaptic neuron
    Panel(v.segments[0].filter(name='v')[0],
          ylabel="Membrane potential (mV)",
          data_labels=[main_pop.label], yticks=True, xlim=(0, runtime)),
    Panel(gsyn_exc.segments[0].filter(name='gsyn_exc')[0],
          ylabel="gsyn excitatory (mV)",
          data_labels=[main_pop.label], yticks=True, xlim=(0, runtime)),
    Panel(gsyn_inh.segments[0].filter(name='gsyn_inh')[0],
          ylabel="gsyn inhibitory (mV)",
          data_labels=[main_pop.label], yticks=True, xlim=(0, runtime)),
    title="Simple synfire chain example",
    annotations="Simulated with {}".format(sim.name())
)
plt.show()
print(figure_filename)

sim.end()
