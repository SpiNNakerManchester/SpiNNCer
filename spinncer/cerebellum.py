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
        self.__number_of_neurons = {k: None for k in POPULATION_ID.keys()}
        self.__nid_offset = {k: None for k in POPULATION_ID.keys()}
        self.__projections = {k: None for k in CONNECTIVITY_MAP.keys()}
        self.__connections = {k: None for k in CONNECTIVITY_MAP.keys()}

        __connectivity = self.__extract_connectivity(connectivity)

        # Construct PyNN neural Populations
        self.__build_populations(np.asarray(__connectivity['positions']))

        # Construct PyNN Projections
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
        # Unique populations ids
        unique_ids = np.unique(positions[:, 1]).astype(int)
        for ui in unique_ids:
            # Get cell name based on the population id
            _cell_name = CELL_NAME_FOR_ID[ui]
            _no_cells = positions[positions[:, 1] == ui, :].shape[0]
            # Store number of neurons for later
            self.__number_of_neurons[_cell_name] = _no_cells
            # Adding the population to the network
            self.__populations[_cell_name] = self.__sim.Population(
                _no_cells, CELL_TYPES[_cell_name], label=_cell_name + " cells")
            # save global neuron ID offset
            # NOTE: np.min(unique_ids) == 1
            if ui == 1:
                self.__nid_offset[_cell_name] = 0
            else:
                self.__nid_offset[_cell_name] = \
                    self.__nid_offset[CELL_NAME_FOR_ID[ui - 1]] + \
                    self.__number_of_neurons[CELL_NAME_FOR_ID[ui - 1]]

    def __build_projections(self, connections):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            projection_reporting(connections)
        # Retrieve the Projection labels
        labels = connections.keys()
        # Loop over each connection in `connections`
        for conn_label in labels:
            # Retrieve saved connectivity and cast to np.ndarray
            # [:, :2] -- the columns are
            # 1. Unique NID for pre-synaptic neuron
            # 2. Unique NID for post-synaptic neuron
            # 3. 3D distance between somas
            conns = np.asarray(connections[conn_label])[:, :2].astype(int)
            no_synapses = conns.shape[0]
            pre_pop = CONNECTIVITY_MAP[conn_label]['pre']
            post_pop = CONNECTIVITY_MAP[conn_label]['post']
            weight = CONNECTIVITY_MAP[conn_label]['weight']
            delay = CONNECTIVITY_MAP[conn_label]['delay']


            # Normalise the source and target neuron IDs
            # Neurons IDs used here are np.arange(0, TOTAL_NUMBER_OF_NEURONS)
            norm_ids = np.asarray([self.__nid_offset[pre_pop],
                                   self.__nid_offset[post_pop]])
            conns -= norm_ids

            # Save the explicit connectivity for later
            stacked_weights = np.asarray([[weight]] * no_synapses)
            stacked_delays = np.asarray([[delay]] * no_synapses)
            self.__connections[conn_label] = np.concatenate([conns,
                                           stacked_weights,
                                            stacked_delays], axis=1)

            # Adding the projection to the network
            self.__projections[conn_label] = self.__sim.Projection(
                self.__populations[pre_pop],  # pre-synaptic population
                self.__populations[post_pop],  # post-synaptic population
                self.__sim.FromListConnector(conns),  # connector
                synapse_type=self.__sim.StaticSynapse(
                    weight=weight,
                    delay=delay),  # synapse type (weights + delays)
                label=conn_label)  # label for connection

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

    @property
    def number_of_neurons(self):
        return self.__number_of_neurons

    @property
    def connections(self):
        return self.__connections
