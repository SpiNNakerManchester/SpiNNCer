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
from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq


class Cerebellum(Circuit):

    def __init__(self, sim, connectivity, stimulus_information, reporting=True,
                 skip_projections=False, weight_scaling=None):
        self.sim = sim
        self.reporting = reporting

        self.populations = {k: None for k in POPULATION_ID.keys()}
        self.number_of_neurons = {k: None for k in POPULATION_ID.keys()}
        self.nid_offset = {k: None for k in POPULATION_ID.keys()}
        self.projections = {k: None for k in CONNECTIVITY_MAP.keys()}
        self.connections = {k: None for k in CONNECTIVITY_MAP.keys()}
        self.stimulus_information = stimulus_information

        self.periodic_stimulus = stimulus_information['periodic_stimulus']
        self.stimulus = None
        self.weight_scaling = weight_scaling or 1.0

        __connectivity = self.__extract_connectivity(connectivity)
        self.cell_positions = np.asarray(__connectivity['positions'])

        # Construct PyNN neural Populations
        self.build_populations(self.cell_positions)

        # Construct PyNN Projections
        if not skip_projections:
            self.build_projections(__connectivity['connections'])

        for pop_name, pop_obj in self.populations.items():
            self.__setattr__(pop_name, pop_obj)

    def build_populations(self, positions):
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
            self.number_of_neurons[_cell_name] = _no_cells
            # save global neuron ID offset
            # NOTE: np.min(unique_ids) == 1
            if ui == 1:
                self.nid_offset[_cell_name] = 0
            else:
                self.nid_offset[_cell_name] = \
                    self.nid_offset[CELL_NAME_FOR_ID[ui - 1]] + \
                    self.number_of_neurons[CELL_NAME_FOR_ID[ui - 1]]
            # Retrieve correct cell parameters for the current cell
            cell_model = CELL_TYPES[_cell_name]
            if _cell_name == "glomerulus":
                cell_param = self.__compute_stimulus(
                    self.stimulus_information, _no_cells)
                additional_params = {'seed': 31415926}
                if self.periodic_stimulus:
                    cell_model = sim.SpikeSourceArray
                    CELL_TYPES[_cell_name] = cell_model
                    additional_params = {}
                self.stimulus = cell_param
            else:
                cell_param = CELL_PARAMS[_cell_name]
                additional_params = {}
                # add E_rev_I to all cells
                if cell_model == self.sim.IF_cond_exp:
                    cell_param['e_rev_I'] = cell_param['v_reset'] - 5.
            # Adding the population to the network
            self.populations[_cell_name] = self.sim.Population(
                _no_cells,
                cellclass=cell_model,
                cellparams=cell_param,
                label=_cell_name + " cells",
                additional_parameters=additional_params)

    def build_projections(self, connections):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            projection_reporting(connections, self.number_of_neurons)
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
                    (CELL_TYPES[post_pop] ==
                     self.sim.extra_models.SpikeSourcePoissonVariable or
                     CELL_TYPES[post_pop] ==
                     self.sim.SpikeSourceArray)):
                # Can't send projections to a spike source
                print("Ignoring connection terminating at", post_pop, "...")
                continue

            # Normalise the source and target neuron IDs
            # Neurons IDs used here are np.arange(0, TOTAL_NUMBER_OF_NEURONS)
            norm_ids = np.asarray([self.nid_offset[pre_pop],
                                   self.nid_offset[post_pop]])
            conns -= norm_ids

            # Save the explicit connectivity for later
            # (after scaling the weights)
            stacked_weights = np.asarray([[np.abs(weight)]] * no_synapses) * \
                self.weight_scaling
            stacked_delays = np.asarray([[delay]] * no_synapses)
            self.connections[conn_label] = np.concatenate(
                [conns, stacked_weights, stacked_delays], axis=1)

            assert (np.max(conns[:, 0]) < self.number_of_neurons[pre_pop]), \
                np.max(conns[:, 0])
            assert (np.max(conns[:, 1]) < self.number_of_neurons[post_pop]), \
                np.max(conns[:, 1])
            # Adding the projection to the network
            receptor_type = "inhibitory" if weight < 0 else "excitatory"
            self.projections[conn_label] = self.sim.Projection(
                self.populations[pre_pop],  # pre-synaptic population
                self.populations[post_pop],  # post-synaptic population
                # connector includes (source, target, weight, delay)
                self.sim.FromListConnector(self.connections[conn_label]),
                receptor_type=receptor_type,  # inh or exc
                label=conn_label)  # label for connection

    def __compute_stimulus(self, stimulus_information, n_inputs,
                           with_positions=True):
        # convert stimulation times to numpy array
        stim_times = np.asarray(stimulus_information['stim_times'])
        f_base = np.asarray(stimulus_information['f_base'])
        f_peak = np.asarray(stimulus_information['f_peak'])
        periodic_stimulus = stimulus_information['periodic_stimulus']

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
        if with_positions:
            radius = self.stimulus_information['stim_radius']
            gloms_pos = self.cell_positions[
                        self.cell_positions[:, 1] == POPULATION_ID['glomerulus'], :]
            # find center of 'glomerular sphere'
            x_c, y_c, z_c = (np.median(gloms_pos[:, 2]),
                             np.median(gloms_pos[:, 3]),
                             np.median(gloms_pos[:, 4]))

            # Find glomeruli falling into the selected volume
            target_gloms_idx = np.sum((gloms_pos[:, 2::] -
                                       np.array([x_c, y_c, z_c])) ** 2,
                                      axis=1).__lt__(radius ** 2)
            target_gloms = np.asarray(gloms_pos[target_gloms_idx, 0]).astype(int)
            # The target_gloms are not normalised (they are global IDs)
            target_gloms -= self.nid_offset['glomerulus']

            # Inverse selection using mask (all other target_gloms are supposed
            # to fire at f_base Hz
            # Thanks to
            # https://stackoverflow.com/questions/25330959/how-to-select-inverse-of-indexes-of-a-numpy-array
            mask = np.ones(n_inputs, np.bool)
            mask[target_gloms] = 0
            # Set the firing rate of other gloms to baseline level
            rates[mask, 1] = f_base
            # Report numbers here
            print("=" * 60)
            print("Number of stimulated Gloms: ", len(target_gloms),
                  "i.e. {:6.2%} the total".format(
                      len(target_gloms) / float(n_inputs)))

        if not periodic_stimulus:
            # VARIABLE RATE POISSON SPIKE SOURCE + INDEPENDENT SPIKE TRAINS
            # if we should only "stimulate" certain cells which I do by setting
            # the rates during stimulation to the base level (f_base)

            # Add a variable-rate poisson spike source to the network
            stimulus_params = {
                'rates': rates,
                'starts': starts,
                'durations': durations,
            }
            return stimulus_params
        else:
            # SPIKE SOURCE ARRAY + PERIODIC STIMULUS
            # prepare spike times for all n_inputs neurons
            spike_times = [[] for _ in range(n_inputs)]
            for i, rate, start, duration in zip(range(n_inputs), rates, starts, durations):
                curr_spikes = []
                for r, s, d in zip(rate, start, duration):
                    if r == f_base:
                        curr_spikes.append(homogeneous_poisson_process(
                            rate=r * pq.Hz,
                            t_start=s * pq.ms,
                            t_stop=(s + d) * pq.ms,
                            as_array=True))
                    else:
                        spike_nums = np.int(np.round((r * d) / 1000.))
                        curr_spikes.append(
                            np.round(np.linspace(s, s + d, spike_nums)))
                spike_times[i] = np.concatenate(curr_spikes)
            return {
                'spike_times': spike_times
            }

    def get_circuit_inputs(self):
        """
        Return a (copy) dictionary of INPUT populations in the Cerebellum
        :return: INPUT populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.populations.items()
                if CELL_IO_STATUS[k] == IO_Status.INPUT}

    def get_circuit_outputs(self):
        """
        Return a (copy) dictionary of OUTPUT populations in the Cerebellum
        :return: OUTPUT populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.populations.items()
                if CELL_IO_STATUS[k] == IO_Status.OUTPUT}

    def get_all_populations(self):
        """
        Return a (copy) dictionary of all populations in the Cerebellum
        :return: populations in the Cerebellar circuit
        :rtype: dict
        """
        return {k: v for k, v in self.populations.items()}

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
        for label, pop in self.populations.items():
            if CELL_TYPES[label] == sim.SpikeSourceArray:
                print("NOT enabling recordings for ", label,
                      "(it's a Spike Source Array)")
                continue
            print("Enabling recordings for ", label, "...")
            pop.record(['spikes'])

    def retrieve_all_recorded_spikes(self):
        """
        Retrieve the recorded spikes for all populations
        :return: spike times for all populations
        :rtype: list or Neo.Block
        """
        all_spikes = {}
        for label, pop in self.populations.items():
            print("Retrieving recordings for ", label, "...")
            if CELL_TYPES[label] == sim.SpikeSourceArray:
                _spikes = []
                for i, _times in enumerate(self.stimulus['spike_times']):
                    for t in _times:
                        _spikes.append(np.asarray([i, t]))
                _spikes = np.asarray(_spikes)
                all_spikes[label] = _spikes
            else:
                all_spikes[label] = pop.get_data(['spikes'])
        return all_spikes

    def retrieve_selective_recordings(self):
        """
        Retrieve the recorded observables for all populations
        :return: spike times for all populations
        :rtype: list or Neo.Block
        """
        gsyn_rec = {}
        for label, pop in self.populations.items():
            if label == "glomerulus":
                print("Skipping selective recording for", label, "...")
                continue
            print("Retrieving recordings for ", label, "...")
            gsyn_rec[label] = {}
            gsyn_rec[label]['gsyn_inh'] = pop.get_data(['gsyn_inh'])
            gsyn_rec[label]['gsyn_exc'] = pop.get_data(['gsyn_exc'])
            gsyn_rec[label]['v'] = pop.get_data(['v'])
        return gsyn_rec

    def selectively_record_all(self, number_of_neurons=None, every=None):
        if bool(number_of_neurons) == bool(every):
            raise ValueError("Specify a number of neurons (sampled randomly) "
                             "to record from xor a linspace of neurons "
                             "(every nth).")
        for label, pop in self.populations.items():
            if label == "glomerulus":
                print("Skipping selective recording for", label, "...")
            else:
                print("Selectively recording gsyn and v for ", label)
                ps = pop.size
                if number_of_neurons:
                    _neuron_choice = np.random.choice(
                        ps, size=min(number_of_neurons, ps), replace=False)
                else:
                    _neuron_choice = np.arange(0, ps, every)
                pop[_neuron_choice].record(['gsyn_inh', 'gsyn_exc', 'v'])

    def retrieve_final_connectivity(self):
        all_connections = {}
        for label, p in self.projections.items():
            if p is None:
                print("Projection", label, "is not implemented!")
                continue
            print("Retrieving connectivity for projection ", label, "...")
            try:
                all_connections[label] = \
                    np.array(p.get(('weight', 'delay'),
                                   format="list")._get_data_items())
            except Exception as e:
                print("Careful! Something happened when retrieving the "
                      "connectivity:", e, "\nRetrying...")
                all_connections[label] = \
                    np.array(p.get(('weight', 'delay'), format="list"))
        return all_connections

    def retrieve_population_names(self):
        return list(self.populations.keys())

    @staticmethod
    def retrieve_cell_params():
        return CELL_PARAMS
