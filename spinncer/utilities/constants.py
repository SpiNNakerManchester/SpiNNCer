from enum import Enum
# import sPyNNaker
try:
    import spynnaker8 as sim
except:
    import pyNN.spynnaker as sim

# Helpful Enum to easily check whether a population should e.g be
# returned from a Circuit object when a user asks for Inputs / Outputs
# or all of the units. An alternative use case: this Enum would control the
# neuron cell type/model. For example, the input population should be
# represented by a new neuron type: Passthrough Neuron that would project with
# both Excitatory and Inhibitory synapses; the output population could be
# represented by a Live Output population / External device.
class IO_Status(Enum):
    OUTPUT = 0
    INPUT = 1
    HIDDEN = 2


# ID is used in hdf5 file defining the positions of individual cells
POPULATION_ID = {
    'golgi': 1,
    'glomerulus': 2,
    'granule': 3,
    'purkinje': 4,
    'basket': 5,
    'stellate': 6,
    'dcn': 7
}

# Cell types
CELL_TYPES = {
    'golgi': sim.IF_cond_exp,
    'glomerulus': sim.IF_cond_exp,
    'granule': sim.IF_cond_exp,
    'purkinje': sim.IF_cond_exp,
    'basket': sim.IF_cond_exp,
    'stellate': sim.IF_cond_exp,
    'dcn': sim.IF_cond_exp
}

# Python 3+ syntax to invert a dictionary
CELL_NAME_FOR_ID = {v: k for k, v in POPULATION_ID.items()}

# values check on 16.10.2019
CELL_PARAMS = {
    'golgi': {'tau_refrac': 2.0,  # ms
              'cm': 0.076,  # nF
              'v_thresh': -55.0,  # mV
              'v_rest': -65.0,  # mV, added 16.10.2019
              'v_reset': -75.0,  # mV
              'tau_m': 21,  # ms, changed from 21.1 on 16.10.2019
              'i_offset': 36.8,  # pA, changed from 36.75 on 16.10.2019
              'tau_syn_E': 0.5,  # ms
              'tau_syn_I': 10.0  # ms
              },
    ' glomerulus': {  # TODO figure out how to make this passthrough

    },  # Glom is special. It's a non-neural mossy fiber terminal (input)
    'granule': {'tau_refrac': 1.5,  # ms
                'cm': 0.003,  # nF
                'v_thresh': -42.0,  # mV
                'v_rest': -74.0,  # mV, added 16.10.2019
                'v_reset': -84.0,  # mV
                'tau_m': 2,  # ms, changed from 20 on 16.10.2019
                'i_offset': 0,  # pA
                'tau_syn_E': 0.5,
                'tau_syn_I': 10.0},
    'purkinje': {'tau_refrac': 0.8,  # ms
                 'cm': 0.62,  # nF, changed from 0.062 on 16.10.2019
                 'v_thresh': -47.0,  # mV
                 'v_rest': -62.0,  # mV, added 16.10.2019
                 'v_reset': -72.0,  # mV
                 'tau_m': 88.0,  # ms, changed from 88.6
                 'i_offset': 600,  # pA, changed from 700
                 'tau_syn_E': 0.5,  # ms
                 'tau_syn_I': 1.6  # ms, changed from 10
                 },
    'basket': {'tau_refrac': 1.6,  # ms, changed from 1.59 on 16.10.2019
               'cm': 0.0146,  # nF
               'v_thresh': -53.0,  # mV
               'v_rest': -68.0,  # mV, added 16.10.2019
               'v_reset': -78.0,  # mV
               'tau_m': 14.6,  # ms
               'i_offset': 15.6,  # pA
               'tau_syn_E': 0.64,  # ms, changed from 0.5 on 16.10.2019
               'tau_syn_I': 2  # ms, changed from 10.0 on 16.10.2019
               },
    'stellate': {'tau_refrac': 1.6,  # ms, changed from 1.59 on 16.10.2019
                 'cm': 0.0146,  # nF
                 'v_thresh': -53.0,  # mV
                 'v_rest': -68.0,  # mV, added 16.10.2019
                 'v_reset': -78.0,  # mV
                 'tau_m': 14.6,  # ms
                 'i_offset': 15.6,  # pA
                 'tau_syn_E': 0.64,  # ms, changed from 0.5 on 16.10.2019
                 'tau_syn_I': 2  # ms, changed from 10.0 on 16.10.2019
                 },
    'dcn': {'tau_refrac': 3.7,  # ms
            'cm': 0.089,  # nF
            'v_thresh': -48.0,  # mV
            'v_rest': -59.0,  # mV, added 16.10.2019
            'v_reset': -69.0,  # mV
            'tau_m': 57,  # ms, changed from 57.1
            'i_offset': 45.75,  # pA
            'tau_syn_E': 7.1,  # ms, changed from 0.5
            'tau_syn_I': 13.6  # ms, changed from 10.0
            }
}

CELL_IO_STATUS = {
    'golgi': IO_Status.HIDDEN,
    'glomerulus': IO_Status.INPUT,
    'granule': IO_Status.HIDDEN,
    'purkinje': IO_Status.HIDDEN,
    'basket': IO_Status.HIDDEN,
    'stellate': IO_Status.HIDDEN,
    'dcn': IO_Status.OUTPUT
}