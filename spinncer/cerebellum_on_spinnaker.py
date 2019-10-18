'''

Script to simulate a model of a cerebellum:
https://www.frontiersin.org/articles/10.3389/fninf.2019.00037/full

This script loads the network connectivity from a hdf5 file and generates
the appropriate PyNN architecture onto which to load this connectivity.
The network is subsequently ran on the SpiNNaker neuromorphic platform, with
simulation results (generally, spikes) saved to a numpy compressed archive.

Analysis of the recorded observables is performed in a different location from
this script.

'''
# imports
import numpy as np, h5py
from spinncer.circuit import Circuit
from spinncer.utilities.constants import *
from colorama import Fore, Style, init as color_init


class Cerebellum(Circuit):

    def __init__(self, sim, connectivity, reporting=True):
        self.__sim = sim
        self.reporting = reporting

        self.__populations = {k: None for k in POPULATION_ID.keys()}
        self.__projections = []

        __connectivity = self.__extract_connectivity(connectivity)

        # Construct PyNN neural Populations
        self.__build_populations(np.asarray(__connectivity['positions']))

        # Construct PyNN Projections
        self.__build_projections(np.asarray(__connectivity['connections']))

        for pop_name, pop_obj in self.__populations.items():
            self.__setattr__(pop_name, pop_obj)

        pass

    def __build_populations(self, positions):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            self.__population_reporting(positions)
        unique_ids = np.unique(positions[:, 1]).astype(int)
        for ui in unique_ids:
            _cell_name = CELL_NAME_FOR_ID[ui]
            _no_cells = positions[positions[:, 1] == ui, :].shape[0]
            # Adding the population to the network
            self.__populations[_cell_name] = sim.Population(
                _no_cells, CELL_TYPES[_cell_name], label=_cell_name + "cells")

    def __build_projections(self, connections):
        """
        Construct PyNN Projections between Populations
        """
        pass

    def get_circuit_inputs(self):
        """
        Return a (copy) dictionary of INPUT populations in the Cerebellum
        :return: INPUT populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k:v for k,v in self.__populations.items()
                if CELL_IO_STATUS[k] == IO_Status.INPUT}


    def get_circuit_outputs(self):
        """
        Return a (copy) dictionary of OUTPUT populations in the Cerebellum
        :return: OUTPUT populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.__populations.items()
                if CELL_IO_STATUS[k] == IO_Status.OUTPUT}

    def get_all_populations(self):
        """
        Return a (copy) dictionary of all populations in the Cerebellum
        :return: populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.__populations.items()}

    def __extract_connectivity(self, connectivity):
        """
        Extract the appropriate information from a dataset of the
        Cerebellum model. It should contain information about the
        number of neurons in each population and the connectivity
        present in the network.
        :param connectivity: file path or Dataset
        :type connectivity: str or Dataset
        :return:
        :rtype: dict
        """
        if isinstance(connectivity, str):
            ext = connectivity.split(".")[-1]
            if ext in ['h5', 'hdf5']:
                return h5py.File(connectivity, 'r')
            else:
                raise NotImplementedError(
                    "The file extension", ext, "is not currently supported.")
        else:
            # We can only hope the data structure is already in an
            # understandable format
            return connectivity

    def __population_reporting(self, positions):
        """
        Helper function reporting on various aspects of the Populations
        :param positions: X Y Z positions of individual cells
        :type positions: np.ndarray
        :return: None
        :rtype: None
        """
        color_init(strip=False)
        unique_ids = np.unique(positions[:, 1]).astype(int)
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
            print("\t{:10} -> {:10} ".format(_cell_name, _no_cells))


if __name__ == "__main__":
    # import sPyNNaker
    try:
        import spynnaker8 as sim
    except:
        import pyNN.spynnaker as sim
    from spinncer_argparser import *
    import pylab as plt

    # Record simulation start time (wall clock)
    start_time = plt.datetime.datetime.now()
    connectivity_filename = 'datasets/scaffold_detailed__158.0x158.0_v3.hdf5'

    # Set up the simulation
    # TODO control some of these using the argparser
    sim.setup(timestep=.1, min_delay=.1, max_delay=10)

    cerebellum_circuit = Cerebellum(sim, connectivity_filename)
    # Test various exposed methods
    populations = cerebellum_circuit.get_all_populations()
    assert (len(populations) == 7)
    input_populations = cerebellum_circuit.get_circuit_inputs()
    assert (len(input_populations) == 1)
    output_populations = cerebellum_circuit.get_circuit_outputs()
    assert (len(output_populations) == 1)

    # Add stimulus to the network

    # Set up recordings

    # Run the simulation

    # Retrieve recordings

    # Compute time taken to reach this point
    end_time = plt.datetime.datetime.now()
    total_time = end_time - start_time
    # Save results

    # Appropriately end the simulation
    sim.end()
    # Report time again
    print(Fore.GREEN + "Total time elapsed -- " + str(total_time) + Style.RESET_ALL)
