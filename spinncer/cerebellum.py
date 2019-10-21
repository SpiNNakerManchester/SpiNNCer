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
import numpy as np
import h5py
from spinncer.circuit import Circuit
from spinncer.utilities.constants import *
from spinncer.utilities.reporting import (population_reporting,
                                          projection_reporting)


class Cerebellum(Circuit):

    def __init__(self, sim, connectivity, reporting=True,
                 skip_projections=False):
        self.__sim = sim
        self.reporting = reporting

        self.__populations = {k: None for k in POPULATION_ID.keys()}
        self.__projections = []

        __connectivity = self.__extract_connectivity(connectivity)

        # Construct PyNN neural Populations
        self.__build_populations(np.asarray(__connectivity['positions']))

        # TODO Construct PyNN Projections
        if not skip_projections:
            self.__build_projections(__connectivity['connections'])

        for pop_name, pop_obj in self.__populations.items():
            self.__setattr__(pop_name, pop_obj)

    def __build_populations(self, positions):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            population_reporting(positions)
        unique_ids = np.unique(positions[:, 1]).astype(int)
        for ui in unique_ids:
            _cell_name = CELL_NAME_FOR_ID[ui]
            _no_cells = positions[positions[:, 1] == ui, :].shape[0]
            # Adding the population to the network
            self.__populations[_cell_name] = self.__sim.Population(
                _no_cells, CELL_TYPES[_cell_name], label=_cell_name + "cells")

    def __build_projections(self, connections):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            projection_reporting(connections)

    def get_circuit_inputs(self):
        """
        Return a (copy) dictionary of INPUT populations in the Cerebellum
        :return: INPUT populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.__populations.items()
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

    def record_all_spikes(self):
        for label, pop in self.__populations.items():
            print("Enabling recordings for ", label, "...")
            pop.record(['spikes'])

    def retrieve_all_recorded_spikes(self, spinnaker_data=True):
        all_spikes = {}
        for label, pop in self.__populations.items():
            print("Retrieving recordings for ", label, "...")
            if spinnaker_data:
                all_spikes[label] = pop.spinnaker_get_data(['spikes'])
            else:
                all_spikes[label] = pop.get_data(['spikes'])
        return all_spikes

    def retrieve_population_names(self):
        return list(self.__populations.keys())

    @staticmethod
    def retrieve_cell_params():
        return CELL_PARAMS
