from colorama import Fore, Style, init as color_init
from spinncer.utilities.constants import (CELL_NAME_FOR_ID,
                                          CELL_IO_STATUS,
                                          IO_Status,
                                          CONNECTIVITY_MAP)
import numpy as np


def population_reporting(positions):
    """
    Helper function reporting on various aspects of the Populations
    :param positions: X Y Z positions of individual cells
    :type positions: np.ndarray
    :return: None
    :rtype: None
    """
    color_init(strip=False)
    unique_ids = np.unique(positions[:, 1]).astype(int)
    total_number_of_neurons = 0
    print("=" * 60)
    print("The file contains information about",
          Fore.GREEN, unique_ids.size, Style.RESET_ALL, "populations")
    print("=" * 60)
    print("Mapping of IDs:")
    print("-" * 60)
    for ui in unique_ids:
        _cell_name = CELL_NAME_FOR_ID[ui]
        _status = CELL_IO_STATUS[_cell_name]
        if _status == IO_Status.INPUT:
            _color = Fore.GREEN
        elif _status == IO_Status.OUTPUT:
            _color = Fore.RED
        else:
            _color = ''
        print("\t{:2d} -> {:10} ".format(ui, _cell_name),
              _color, "[{:16}]".format(_status), Style.RESET_ALL)
    print("=" * 60)
    print("Number of neurons in each population")
    print("-" * 60)
    for ui in unique_ids:
        _cell_name = CELL_NAME_FOR_ID[ui]
        _no_cells = positions[positions[:, 1] == ui, :].shape[0]
        total_number_of_neurons += _no_cells
        print("\t{:10} -> {:10} neurons".format(_cell_name, _no_cells))
    print(Fore.GREEN, "\t{:10} -> {:10} neurons".format(
        "TOTAL", total_number_of_neurons), Style.RESET_ALL)


def projection_reporting(connections):
    """
    Helper function reporting on various aspects of the Projections
    :param connections: pre-id, post-id, xxxx
    :type connections: HDF5 group
    :return: None
    :rtype: None
    """
    color_init(strip=False)
    connection_keys = connections.keys()
    number_of_afferents = {k: 0 for k in CELL_IO_STATUS.keys()}
    total_number_of_synapses = 0
    print("=" * 60)
    print("The file contains information about",
          Fore.GREEN, len(connection_keys), Style.RESET_ALL, "projections")
    print("=" * 60)
    print("Number of synapses per projection:")
    print("-" * 60)
    for key in connection_keys:
        # report number of synapses
        conns = np.asarray(connections[key])
        number_of_syn = conns.shape[0]
        total_number_of_synapses += number_of_syn
        print("\t{:10} -> {:10} synapses".format(key, number_of_syn))
        # compute fan in of post population
        number_of_afferents[CONNECTIVITY_MAP[key]['post']] += number_of_syn
    print(Fore.GREEN, "\t{:10} -> {:10} synapses".format(
        "TOTAL", total_number_of_synapses), Style.RESET_ALL)
    print("=" * 60)
    print("Number of incoming connections per population:")
    print("-" * 60)
    for pop, fan_in in number_of_afferents.items():
        print("\t{:10} -> {:10} incoming synapses".format(pop, fan_in))
    print("=" * 60)
