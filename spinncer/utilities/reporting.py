from colorama import Fore, Style, init as color_init
from spinncer.utilities.constants import (CELL_NAME_FOR_ID,
                                          CELL_IO_STATUS,
                                          IO_Status,
                                          CONNECTIVITY_MAP)
import numpy as np


def population_reporting(positions, number_of_neurons):
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
    if unique_ids.size == len(CELL_NAME_FOR_ID.keys()):
        print("=" * 80)
        print("The file contains information about",
              Fore.GREEN, unique_ids.size, Style.RESET_ALL, "populations")
        print("=" * 80)
        print("Mapping of IDs:")
        print("-" * 80)
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
        print("=" * 80)
        print("Number of neurons in each population")
        print("-" * 80)
        for ui in unique_ids:
            _cell_name = CELL_NAME_FOR_ID[ui]
            _no_cells = positions[positions[:, 1] == ui, :].shape[0]
            total_number_of_neurons += _no_cells
            print("\t{:15} -> {:10} neurons".format(_cell_name, _no_cells))
        print(Fore.GREEN, "\t{:10} -> {:10} neurons".format(
            "TOTAL", total_number_of_neurons), Style.RESET_ALL)
    else:
        print("=" * 80)
        print("Number of neurons in each population")
        print("-" * 80)
        for _cell_name, _no_cells in number_of_neurons.items():
            total_number_of_neurons += _no_cells
            print("\t{:15} -> {:10} neurons".format(_cell_name, _no_cells))
        print(Fore.GREEN, "\t{:10} -> {:10} neurons".format(
            "TOTAL", total_number_of_neurons), Style.RESET_ALL)
        print("=" * 80)


def projection_reporting(connections, number_of_neurons, conn_params):
    """
    Helper function reporting on various aspects of the Projections
    :param connections: pre-id, post-id, xxxx
    :type connections: HDF5 group
    :return: None
    :rtype: None
    """
    color_init(strip=False)
    connection_keys = conn_params.keys()
    number_of_afferents = {k: 0 for k in number_of_neurons.keys()}
    total_number_of_synapses = 0
    print("=" * 80)
    print("The file contains information about",
          Fore.GREEN, len(connection_keys), Style.RESET_ALL, "projections")
    print("=" * 80)
    print("Number of synapses per projection:")
    print("-" * 80)
    for key in connection_keys:
        # report number of synapses
        conns = np.asarray(connections[key])
        weight = conn_params[key]['weight']
        if weight >= 0:
            coloured_syn_type = Fore.GREEN + "[exc]" + Style.RESET_ALL
        else:
            coloured_syn_type = Fore.RED + "[inh]" + Style.RESET_ALL
        number_of_syn = conns.shape[0]
        total_number_of_synapses += number_of_syn
        print("\t{:27} -> {:10} synapses {:5}".format(key,
                                                      number_of_syn,
                                                      coloured_syn_type))
        # compute fan in of post population
        number_of_afferents[conn_params[key]['post']] += number_of_syn
    print(Fore.GREEN, "\t{:10} -> {:10} synapses".format(
        "TOTAL", total_number_of_synapses), Style.RESET_ALL)
    print("=" * 80)
    print("Number of incoming connections per population:")
    print("-" * 80)
    for pop, fan_in in number_of_afferents.items():
        print("\t{:15} -> {:10} incoming synapses".format(pop, fan_in))
    print("=" * 80)
    print("Normalised number of incoming connections per population:")
    print("-" * 80)
    for pop, fan_in in number_of_afferents.items():
        print("\t{:15} -> {:>8.2f} incoming synapses".format(
            pop, fan_in / number_of_neurons[pop]))
    print("=" * 80)
