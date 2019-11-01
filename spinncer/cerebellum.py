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

    def __init__(self, sim, connectivity, stimulus_information, reporting=True,
                 skip_projections=False):
        self.__sim = sim
        self.reporting = reporting

        self.__populations = {k: None for k in POPULATION_ID.keys()}
        self.__number_of_neurons = {k: None for k in POPULATION_ID.keys()}
        self.__nid_offset = {k: None for k in POPULATION_ID.keys()}
        self.__projections = {k: None for k in CONNECTIVITY_MAP.keys()}
        self.__connections = {k: None for k in CONNECTIVITY_MAP.keys()}
        self.__stimulus_information = stimulus_information

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
            # Retrieve correct cell parameters for the current cell
            if _cell_name == "glomerulus":
                cell_param = self.compute_stimulus(
                    self.__stimulus_information, _no_cells)
                additional_params = {'seed': 31415926}
            else:
                cell_param = CELL_PARAMS[_cell_name]
                additional_params = {}
            # Adding the population to the network
            self.__populations[_cell_name] = self.__sim.Population(
                _no_cells,
                cellclass=CELL_TYPES[_cell_name],
                cellparams=cell_param,
                label=_cell_name + " cells",
                additional_parameters=additional_params)
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
            projection_reporting(connections, self.__number_of_neurons)
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
            print("Creating projection from {:10}".format(pre_pop),
                  "to {:10}".format(post_pop),
                  "with a weight of {: 2.6f}".format(weight),
                  "uS and a delay of", delay, "ms")
            if (post_pop == "glomerulus" and
                    CELL_TYPES[post_pop] ==
                    self.__sim.extra_models.SpikeSourcePoissonVariable):
                # Can't send projections to a spike source
                continue

            # Normalise the source and target neuron IDs
            # Neurons IDs used here are np.arange(0, TOTAL_NUMBER_OF_NEURONS)
            norm_ids = np.asarray([self.__nid_offset[pre_pop],
                                   self.__nid_offset[post_pop]])
            conns -= norm_ids

            # Save the explicit connectivity for later
            stacked_weights = np.asarray([[np.abs(weight)]] * no_synapses)
            stacked_delays = np.asarray([[delay]] * no_synapses)
            self.__connections[conn_label] = np.concatenate(
                [conns, stacked_weights, stacked_delays], axis=1)

            assert (np.max(conns[:, 0]) < self.number_of_neurons[pre_pop]), \
                np.max(conns[:, 0])
            assert (np.max(conns[:, 1]) < self.number_of_neurons[post_pop]), \
                np.max(conns[:, 1])
            # Adding the projection to the network
            receptor_type = "inhibitory" if weight < 0 else "excitatory"
            self.__projections[conn_label] = self.__sim.Projection(
                self.__populations[pre_pop],  # pre-synaptic population
                self.__populations[post_pop],  # post-synaptic population
                # connector includes (source, target, weight, delay)
                self.__sim.FromListConnector(self.__connections[conn_label]),
                receptor_type=receptor_type,  # inh or exc
                label=conn_label)  # label for connection

    @staticmethod
    def compute_stimulus(stimulus_information, n_inputs):
        # convert stimulation times to numpy array
        stim_times = np.asarray(stimulus_information['stim_times'])
        f_base = np.asarray(stimulus_information['f_base'])
        f_peak = np.asarray(stimulus_information['f_peak'])
        # compute number of rate changes required
        number_of_slots = int(stim_times.size)
        # compute the time at which individual rates take effect
        stim_start = np.cumsum(np.append([0], stim_times))[:number_of_slots]
        starts = np.ones((n_inputs, number_of_slots)) * stim_start
        # compute the duration (in ms) for which individual rates are active
        durations = np.ones((n_inputs, number_of_slots)) * stim_times
        # compute the individual rates (in Hz) for each slot
        rates = np.ones((n_inputs, number_of_slots)) * \
                np.asarray([f_base, f_base + f_peak, f_base])
        # Add a variable-rate poisson spike source to the network
        stimulus_params = {
            'rates': rates,
            'starts': starts,
            'durations': durations,
        }
        return stimulus_params

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
        """
        Retrieve the recorded spikes for all populations
        :param spinnaker_data: if True will return spikes in a 2D list where
        the first column is the neuron id and the second column is the spike
        time of that cell; else spikes will be returned inside a Neo object
        :type spinnaker_data: bool
        :return: spike times for all populations
        :rtype: list or Neo.Block
        """
        all_spikes = {}
        for label, pop in self.__populations.items():
            print("Retrieving recordings for ", label, "...")
            if spinnaker_data:
                all_spikes[label] = pop.spinnaker_get_data(['spikes'])
            else:
                all_spikes[label] = pop.get_data(['spikes'])
        return all_spikes

    def retrieve_final_connectivity(self):
        all_connections = {}
        for label, p in self.__projections.items():
            if p is None:
                print("Projection", label, "is not implemented!")
                continue
            print("Retrieving connectivity for projection ", label, "...")
            all_connections[label] = \
                np.array(p.get(('weight', 'delay'), format="list")._get_data_items())
        return all_connections

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
        """
        This property will return the INITIAL weights
        :return: INITIAL weights
        :rtype: dict
        """
        return self.__connections
