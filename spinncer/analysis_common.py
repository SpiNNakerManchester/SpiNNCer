import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm as cm_mlib
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy
from matplotlib import animation, rc, colors
from brian2.units import *
import matplotlib as mlib
from scipy import stats
from pprint import pprint as pp
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
import traceback
import os
import copy
import neo
from datetime import datetime
import warnings
import ntpath
from spinncer.utilities.constants import CONNECTIVITY_MAP, CELL_PARAMS
from spinncer.utilities.neo_converter import convert_spikes
from colorama import Fore, Style, init as color_init
import pandas as pd
import string
from matplotlib.ticker import MultipleLocator

mlib.use('Agg')
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ensure we use viridis as the default cmap
plt.viridis()

# ensure we use the same rc parameters for all matplotlib outputs
mlib.rcParams.update({'font.size': 24})
mlib.rcParams.update({'errorbar.capsize': 5})
mlib.rcParams.update({'figure.autolayout': True})
viridis_cmap = mlib.cm.get_cmap('viridis')

PREFFERED_ORDER = ['mossy_fibers',
                   'glomerulus',
                   'granule', 'golgi',
                   'stellate', 'basket',
                   'purkinje', 'dcn']


def color_for_index(index, size, cmap=viridis_cmap):
    return cmap(index / (size + 1))


def write_sep():
    print("=" * 80)


def write_line():
    print("-" * 80)


def write_header(msg):
    write_sep()
    print(msg)
    write_line()


def write_short_msg(msg, value):
    print("{:40}:{:39}".format(msg, str(value)))


def write_value(msg, value):
    print("{:60}:{:19}".format(msg, str(value)))


def get_plot_order(for_keys):
    # Compute plot order
    plot_order = []
    # only focus on keys for pops that have spikes
    key_duplicate = list(for_keys)
    key_duplicate.sort()
    for pref in PREFFERED_ORDER:
        for i, key in enumerate(key_duplicate):
            if pref in key:
                plot_order.append(key)
                key_duplicate.pop(i)

    # add remaining keys to the end of plot order
    plot_order += key_duplicate
    return plot_order


COMMON_DISPLAY_NAMES = {
    'f_peak': "$f_{peak}$ (Hz)",
    'spinnaker': "SpiNNaker",
    'nest': "NEST",
    'stim_radius': "Stimulation radius ($\mu m$)",
    'glomerulus cells': "Glom",
    'glomerulus': "Glom",
    'granule cells': "GrC",
    'granule': "GrC",
    'dcn cells': "DCNC",
    'dcn': "DCNC",
    'golgi cells': "GoC",
    'golgi': "GoC",
    'purkinje cells': "PC",
    'purkinje': "PC",
    'stellate cells': "SC",
    'stellate': "SC",
    'basket cells': "BC",
    'basket': "BC",
    'max_spikes_in_a_tick': "Max # of spikes in a time step",
    'send_multicast_packets': "Multicast packets sent",
    'router_provenance': "Router"
}

def capitalise(name):
    return string.capwords(
        " ".join(
            str.split(name, "_")
        ))

def use_display_name(name):
    name = name.lower()
    return COMMON_DISPLAY_NAMES[name] \
        if name in COMMON_DISPLAY_NAMES.keys() \
        else capitalise(name)


def save_figure(plt, name, extensions=(".png", ), **kwargs):
    for ext in extensions:
        write_short_msg("Plotting", name+ext)
        plt.savefig(name + ext, **kwargs)
