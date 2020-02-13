from enum import Enum


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
_model = "IF_cond_exp"
CELL_TYPES = {
    'golgi': _model,
    'glomerulus': "SpikeSourcePoissonVariable",
    'granule': _model,
    'purkinje': _model,
    'basket': _model,
    'stellate': _model,
    'dcn': _model
}

# Python 3+ syntax to invert a dictionary
CELL_NAME_FOR_ID = {v: k for k, v in POPULATION_ID.items()}

# values check on 16.10.2019
# Checked against:
# https://www.frontiersin.org/files/Articles/444802/fninf-13-00037-HTML-r1/image_m/fninf-13-00037-t002.jpg
CELL_PARAMS = {
    'golgi': {'tau_refrac': 2.0,  # ms
              'cm': 0.076,  # nF
              'v_thresh': -55.0,  # mV
              'v_rest': -65.0,  # mV, added 16.10.2019
              'v_reset': -75.0,  # mV
              'tau_m': 21,  # ms, changed from 21.1 on 16.10.2019
              'i_offset': 36.8e-3,  # pA, changed from 36.75 on 16.10.2019
              'tau_syn_E': 0.5,  # ms
              'tau_syn_I': 10.0  # ms
              },
    'glomerulus': {  # TODO figure out how to make this passthrough

    },  # Glom is special. It's a non-neural mossy fiber terminal (input)
    'granule': {'tau_refrac': 1.5,  # ms
                'cm': 0.003,  # nF
                'v_thresh': -42.0,  # mV
                'v_rest': -74.0,  # mV, added 16.10.2019
                'v_reset': -84.0,  # mV
                'tau_m': 2,  # ms, changed from 20 on 16.10.2019
                'i_offset': 0e-3,  # pA
                'tau_syn_E': 0.5,
                'tau_syn_I': 10.0},
    'purkinje': {'tau_refrac': 0.8,  # ms
                 'cm': 0.62,  # nF, changed from 0.062 on 16.10.2019
                 'v_thresh': -47.0,  # mV
                 'v_rest': -62.0,  # mV, added 16.10.2019
                 'v_reset': -72.0,  # mV
                 'tau_m': 88.0,  # ms, changed from 88.6
                 'i_offset': 600e-3,  # pA, changed from 700
                 'tau_syn_E': 0.5,  # ms
                 'tau_syn_I': 1.6  # ms, changed from 10
                 },
    'basket': {'tau_refrac': 1.6,  # ms, changed from 1.59 on 16.10.2019
               'cm': 0.0146,  # nF
               'v_thresh': -53.0,  # mV
               'v_rest': -68.0,  # mV, added 16.10.2019
               'v_reset': -78.0,  # mV
               'tau_m': 14.6,  # ms
               'i_offset': 15.6e-3,  # pA
               'tau_syn_E': 0.64,  # ms, changed from 0.5 on 16.10.2019
               'tau_syn_I': 2  # ms, changed from 10.0 on 16.10.2019
               },
    'stellate': {'tau_refrac': 1.6,  # ms, changed from 1.59 on 16.10.2019
                 'cm': 0.0146,  # nF
                 'v_thresh': -53.0,  # mV
                 'v_rest': -68.0,  # mV, added 16.10.2019
                 'v_reset': -78.0,  # mV
                 'tau_m': 14.6,  # ms
                 'i_offset': 15.6e-3,  # pA
                 'tau_syn_E': 0.64,  # ms, changed from 0.5 on 16.10.2019
                 'tau_syn_I': 2  # ms, changed from 10.0 on 16.10.2019
                 },
    'dcn': {'tau_refrac': 3.7,  # ms
            'cm': 0.089,  # nF
            'v_thresh': -48.0,  # mV
            'v_rest': -59.0,  # mV, added 16.10.2019
            'v_reset': -69.0,  # mV
            'tau_m': 57,  # ms, changed from 57.1
            'i_offset': 55.8e-3,  # pA
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

# PyNN units for synapses:
# http://neuralensemble.org/docs/PyNN/reference/plasticitymodels.html#weight-dependence-components
# suggests that conductance-based synapses have weight = uS
# NB paper weights are defined using nS:
# https://www.frontiersin.org/files/Articles/444802/fninf-13-00037-HTML-r1/image_m/fninf-13-00037-t003.jpg
# delay = ms
CONNECTIVITY_MAP = {
    'aa_goc': {
        'pre': 'granule',
        'post': 'golgi',
        'weight': 20.0e-03,  # uS
        'delay': 2.0,  # ms
    },
    'aa_pc': {
        'pre': 'granule',
        'post': 'purkinje',
        'weight': 75.0e-03,  # uS
        'delay': 2.0,  # ms
    },
    'bc_pc': {
        'pre': 'basket',
        'post': 'purkinje',
        'weight': -9.0e-03,  # uS
        'delay': 4.0,  # ms
    },
    'gj_bc': {
        'pre': 'basket',
        'post': 'basket',
        'weight': -2.5e-03,  # uS
        'delay': 4.0,  # ms
    },
    'gj_goc': {
        'pre': 'golgi',
        'post': 'golgi',
        'weight': -8.0e-03,  # uS
        'delay': 1.0,  # ms
    },
    'gj_sc': {
        'pre': 'stellate',
        'post': 'stellate',
        'weight': -2.0e-03,  # uS
        'delay': 1.0,  # ms
    },
    'glom_dcn': {
        'pre': 'glomerulus',
        'post': 'dcn',
        'weight': 0.006e-03,  # uS
        'delay': 4.0,  # ms
    },
    'glom_goc': {
        'pre': 'glomerulus',
        'post': 'golgi',
        'weight': 2.0e-03,  # uS
        'delay': 4.0,  # ms
    },
    'glom_grc': {
        'pre': 'glomerulus',
        'post': 'granule',
        'weight': 9.0e-03,  # uS
        'delay': 4.0,  # ms
    },
    'goc_grc': {
        'pre': 'golgi',
        'post': 'granule',
        'weight': -5.0e-03,  # uS
        'delay': 2.0,  # ms
    },
    'pc_dcn': {
        'pre': 'purkinje',
        'post': 'dcn',
        'weight': -0.03e-03,  # uS
        'delay': 4.0,  # ms
    },
    'pf_bc': {
        'pre': 'granule',
        'post': 'basket',
        'weight': 0.2e-03,  # uS
        'delay': 5.0,  # ms
    },
    'pf_goc': {
        'pre': 'granule',
        'post': 'golgi',
        'weight': 0.4e-03,  # uS
        'delay': 5.0,  # ms
    },
    'pf_pc': {
        'pre': 'granule',
        'post': 'purkinje',
        'weight': 0.02e-03,  # uS
        'delay': 5.0,  # ms
    },
    'pf_sc': {
        'pre': 'granule',
        'post': 'stellate',
        'weight': 0.2e-03,  # uS
        'delay': 5.0,  # ms
    },
    'sc_pc': {
        'pre': 'stellate',
        'post': 'purkinje',
        'weight': -8.5e-03,  # uS
        'delay': 5.0,  # ms
    },
}

# create a map from pyNest param names to PyNN param names + unit conversions
# https://nest-test.readthedocs.io/en/pynest_mock/models/neurons.html#_CPPv4N4nest14iaf_cond_alphaE
# http://neuralensemble.org/docs/PyNN/standardmodels.html?highlight=if_cond_alpha
PYNEST_TO_PYNN_CONVERSION = {
    "V_m": (None, 0),  # mV
    "E_L": ("v_rest", 1),  # mV
    "C_m": ("cm", 1e-3),  # conversion from pF (pyNest) to nF (PyNN)
    "t_ref": ("tau_refract", 1),  # ms
    "V_th": ("v_thresh", 1),  # mV
    "V_reset": ("v_reset", 1),  # mV
    "E_ex": ("e_rev_E", 1),  # mV
    "E_in": ("e_rev_I", 1),  # mV
    "g_L": (None, 1e-3),  # conversion from nS (pyNest) to uS (PyNN)
    "tau_syn_ex": ("tau_syn_E", 1),  # ms
    "tau_syn_in": ("tau_syn_I", 1),  # ms
    "I_e": ("i_offset", 1e-3),  # conversion from pA (pyNest) to nA (PyNN)
}
