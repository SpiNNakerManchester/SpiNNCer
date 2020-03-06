"""

Script to simulate a model of a cerebellum:
https://www.frontiersin.org/articles/10.3389/fninf.2019.00037/full

This script loads the network connectivity from a hdf5 file and generates
the appropriate PyNN architecture onto which to load this connectivity.
The network is subsequently ran on the SpiNNaker neuromorphic platform, with
simulation results (generally, spikes) saved to a numpy compressed archive.

Analysis of the recorded observables is performed in a different location from
this script.

"""
# imports
import numpy as np
import h5py
from spinncer.circuit import Circuit
from spinncer.utilities.constants import *
from spinncer.utilities.utils import flatten_dict, create_poisson_spikes
from spinncer.utilities.reporting import (population_reporting,
                                          projection_reporting)
from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq
import sys


class Cerebellum(Circuit):

    def __init__(self, sim, connectivity, stimulus_information, reporting=True,
                 params=None, skip_projections=False,
                 weight_scaling=None, save_conversion_file=False,
                 neuron_model="IF_cond_exp", force_number_of_neurons=None,
                 input_spikes=None, rb_left_shifts=None):
        """
        Cerebellum Circuit
        """
        # Reference to the PyNN simulation that this object (Cerebellum) is
        # building itself into
        self.sim = sim
        # Flag that controls whether reports are printed as the  network is
        # being generated
        self.reporting = reporting

        # flag to be used when building the connectivity from hdf5 files
        self.new_scaffold = False
        # Attempt to use passed in cell and connectivity params. If none,
        # use defaults specified in constants but report DECISION
        print("=" * 80)

        connectivity_data = self.__extract_connectivity(connectivity)
        if 'connections' in connectivity_data.keys():
            self.cell_positions = np.asarray(connectivity_data['positions'])
            connections = connectivity_data['connections']
            self.new_scaffold = False
            print("[{:10}]: successfully retrieved connectivity for OLD style"
                  "of scaffold!".format("INFO"))
        else:
            self.cell_positions = np.asarray(connectivity_data['cells']['positions'])
            connections = connectivity_data['cells']['connections']
            self.new_scaffold = True

        if params:
            # Convert from pyNest to PyNN
            # https://nest-test.readthedocs.io/en/pynest_mock/models/neurons.html#_CPPv4N4nest14iaf_cond_alphaE
            # use conversion map from pyNest param names to PyNN param names + units
            # ignore V_m in pyNest
            print("Cell params extracted and converted from Json. ")
            self.cell_params = {}
            sim_params = params['simulations']['FCN_2019']
            for cell_name, param_sets in sim_params['cell_models'].items():
                self.cell_params[cell_name] = {}
                # WARNING this assumes that the different models included in
                # the nested dictionary do NOT have parameter names in common
                # otherwise, "the last one" is the one used here
                if "eglif_cond_alpha_multisyn" in param_sets.keys():
                    del param_sets["eglif_cond_alpha_multisyn"]
                param_dict = flatten_dict(param_sets)
                cp = self.cell_params[cell_name]
                if ('neuron_model' in param_dict.keys() and
                        param_dict['neuron_model'] == "parrot_neuron"):
                    # Mark the current population as input in a distinctive
                    # fashion
                    self.cell_params[cell_name] = None
                    continue
                for p, v in param_dict.items():
                    if (p in PYNEST_TO_PYNN_CONVERSION.keys() and
                            PYNEST_TO_PYNN_CONVERSION[p][0] is not None):
                        cp[PYNEST_TO_PYNN_CONVERSION[p][0]] = v * PYNEST_TO_PYNN_CONVERSION[p][1]
                    if p == "g_L":
                        # convert g_L to tau_m using the equation tau_m = R * C
                        cp['tau_m'] = param_dict['C_m'] / v

            print("Connection params extracted and converted from Json.")
            self.conn_params = {}
            for conn_name, param_sets in sim_params['connection_models'].items():
                # Adjust
                self.conn_params[conn_name] = {}
                adjust_for_units = 1e-3 if sim_params['simulator'] == "nest" else 1
                self.conn_params[conn_name]['weight'] = param_sets['connection']['weight'] * adjust_for_units
                self.conn_params[conn_name]['delay'] = param_sets['connection']['delay']

                # Read connectivity from hdf5 files first rather than from the
                # JSON
                for conn_type in connections[conn_name].attrs['connection_types']:
                    if "from_cell_types" in connections[conn_name].attrs:
                        self.conn_params[conn_name]['pre'] = \
                            connections[conn_name].attrs["from_cell_types"][0]
                    else:
                        self.conn_params[conn_name]['pre'] = \
                            params['connection_types'][conn_type]['from_cell_types'][0]['type']

                    self.conn_params[conn_name]['post'] = \
                        params['connection_types'][conn_type]['to_cell_types'][0]['type']
                    print("[{:8}]: retrieved pre and post for {:15}".format(
                        "INFO", conn_name))
                if self.conn_params[conn_name]['pre'] == "mossy_fibers":
                    self.cell_params["mossy_fibers"] = None
        else:
            print("Cell params not specified. Using defaults...")
            self.cell_params = CELL_PARAMS
            print("Connection params not specified. Using defaults...")
            self.conn_params = CONNECTIVITY_MAP

        print("=" * 80)

        self.populations = {k: None for k in self.cell_params.keys()}
        self.number_of_neurons = {k: None for k in self.cell_params.keys()}
        self.nid_offset = {k: None for k in self.cell_params.keys()}
        self.rb_shifts = {k: None for k in self.cell_params.keys()}
        for k in self.rb_shifts.keys():
            self.rb_shifts[k] = rb_left_shifts
        # Save the neuron model to be used by spiking neurons in the network
        self.neuron_models = {k: str.lower(neuron_model) for k in self.cell_params.keys()}
        # Hard-code glomerulus
        # TODO do better
        if 'glomerulus' in self.neuron_models.keys():
            self.neuron_models['glomerulus'] = "spikesourcearray"
        # Hard-code mossy fibers
        # TODO do better
        if 'mossy_fibers' in self.neuron_models.keys():
            self.neuron_models['mossy_fibers'] = "spikesourcearray"
            # self.number_of_neurons['mossy_fibers'] = 1

        # Force number of neurons
        self.force_number_of_neurons = force_number_of_neurons
        if isinstance(force_number_of_neurons, dict):
            for k, v in force_number_of_neurons.items():
                self.number_of_neurons[k] = v
        elif force_number_of_neurons is not None:
            for k, _ in self.number_of_neurons.items():
                self.number_of_neurons[k] = force_number_of_neurons

        self.projections = {k: None for k in self.conn_params.keys()}
        self.connections = {k: None for k in self.conn_params.keys()}

        self.stimulus_information = stimulus_information

        self.periodic_stimulus = stimulus_information['periodic_stimulus']
        self.weight_scaling = weight_scaling or 1.0
        if self.new_scaffold:
            # Try to read in nid_offsets too
            # populate offsets here too based on ['cells']['type_maps']
            # populate number of neurons here too based on ['cells']['type_maps']
            for k, v in self.nid_offset.items():
                # get the minimum GID to use as an offset
                assert (v is None), "nid_offset key:{} value:{}".format(k, v)
                # Check if cell_name (k) is an entity
                if k in connectivity_data['entities'].keys():
                    actual_ids = np.asarray(connectivity_data['entities'][k])
                else:
                    tm = np.asarray(connectivity_data['cells']['type_maps'][k + "_map"])
                    actual_ids = self.cell_positions[tm][:, 0]
                if actual_ids.size > 0:
                    self.nid_offset[k] = \
                        np.min(actual_ids).astype(int)
                else:
                    # Check if this entity is in the entitites dataset thing
                    self.nid_offset[k] = 0
                # get the number of neurons
                # if k != "mossy_fibers":
                if not force_number_of_neurons:
                    assert (self.number_of_neurons[k] is None), \
                        "number_of_neurons key:{} value:{}".format(
                            k, self.number_of_neurons[k])
                    self.number_of_neurons[k] = \
                        actual_ids.size
            print("[{:10}]: successfully retrieved connectivity for NEW style"
                  "of scaffold!".format("INFO"))
        else:  # old scaffold
            # Compute the number of neurons in each population
            self.__compute_number_of_neurons()

        self._raw_connectivity_info = connections
        if not skip_projections and force_number_of_neurons is None:
            self.__normalise_connections()

        if input_spikes:
            # use passed in spikes
            self.stimulus = input_spikes
        else:
            # generate necessary spikes
            self.stimulus = self.__compute_stimulus()

        # Construct PyNN neural Populations
        self.build_populations(self.cell_positions)
        # Construct PyNN Projections
        if not skip_projections and force_number_of_neurons is None:
            self.build_projections()
        else:
            print("NOT GENERATING ANY CONNECTIVITY FOR THIS MODEL")
            print("skip_projections", skip_projections)
            print("force_number_of_neurons", force_number_of_neurons)

        if save_conversion_file:
            np.savez_compressed("conversion_constants",
                                nid_offsets=self.nid_offset,
                                connectivity=self.connections,
                                all_neurons=self.number_of_neurons,
                                connectivity_file=connectivity)
            print("Saved conversion_constants.npz. Exiting...")
            sys.exit(0)

        for pop_name, pop_obj in self.populations.items():
            self.__setattr__(pop_name, pop_obj)

    def __compute_number_of_neurons(self):
        # compute number of neurons here if not done previously
        unique_ids = np.unique(self.cell_positions[:, 1]).astype(int)
        if unique_ids.size == len(CELL_NAME_FOR_ID.keys()):
            for ui in unique_ids:
                # Get cell name based on the population id
                cell_name = CELL_NAME_FOR_ID[ui]
                no_cells = self.cell_positions[self.cell_positions[:, 1] == ui, :].shape[0]
                # Store number of neurons for later
                self.number_of_neurons[cell_name] = no_cells
                # save global neuron ID offset
                # NOTE: np.min(unique_ids) == 1
                if ui == 1:
                    self.nid_offset[cell_name] = 0
                else:
                    self.nid_offset[cell_name] = \
                        self.nid_offset[CELL_NAME_FOR_ID[ui - 1]] + \
                        self.number_of_neurons[CELL_NAME_FOR_ID[ui - 1]]

    def build_populations(self, positions):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            population_reporting(positions, self.number_of_neurons)

        for cell_name, cell_param in self.cell_params.items():
            # Retrieve correct cell parameters for the current cell
            cell_model = self.neuron_models[cell_name]
            no_cells = self.number_of_neurons[cell_name]
            # skip the creation of glom or mossy fibers
            ins = ["glomerulus", "mossy_fibers"]
            if cell_name in ins and self.force_number_of_neurons:
                print("SKIPPING THE CREATION OF", cell_name,
                      "BECAUSE FORCING # OF NEURONS")
                continue
            if cell_name in ["glomerulus", "mossy_fibers"]:
                cell_param = self.stimulus[cell_name]
                cell_model = self.sim.SpikeSourceArray
                CELL_TYPES[cell_name] = cell_model
                additional_params = {}
            else:
                # else for all other cells
                additional_params = {"rb_left_shifts":
                                     self.rb_shifts[cell_name]}
                # add E_rev_I to all cells
                capp_rev = -90.
                if cell_model == "if_cond_exp":
                    cell_model = self.sim.IF_cond_exp
                    cell_param['e_rev_I'] = capp_rev  # mV
                    cell_param['e_rev_E'] = 0.  # mV
                elif cell_model == "if_curr_exp":
                    cell_model = self.sim.IF_curr_exp
                elif cell_model == "if_curr_alpha":
                    cell_model = self.sim.IF_curr_alpha
                elif cell_model == "if_cond_alpha":
                    cell_model = self.sim.IF_cond_alpha
                    cell_param['e_rev_I'] = capp_rev  # mV
                    cell_param['e_rev_E'] = 0.  # mV
            # Adding the population to the network
            try:
                self.populations[cell_name] = self.sim.Population(
                    no_cells,
                    cellclass=cell_model,
                    cellparams=cell_param,
                    label=cell_name + " cells",
                    additional_parameters=additional_params)
            except TypeError as te:
                self.populations[cell_name] = self.sim.Population(
                    no_cells,
                    cellclass=cell_model,
                    cellparams=cell_param,
                    label=cell_name + " cells")

    def __normalise_connections(self):
        connections = self._raw_connectivity_info
        # Retrieve the Projection labels
        labels = connections.keys()
        # Loop over each connection in `connections`
        for conn_label in labels:
            if conn_label not in self.conn_params.keys():
                print("[{:10}]: CONNECTION UNREGONISED -> "
                      "{:25}".format("ERROR", conn_label))
                continue
            conns = np.asarray(connections[conn_label])[:, :2].astype(int)
            no_synapses = conns.shape[0]
            pre_pop = self.conn_params[conn_label]['pre']
            post_pop = self.conn_params[conn_label]['post']
            weight = self.conn_params[conn_label]['weight']
            delay = self.conn_params[conn_label]['delay']
            if (post_pop == "glomerulus" and pre_pop != "mossy_fibers"):
                # Can't send projections to a spike source
                print("Ignoring connection {:25} "
                      "terminating at".format(conn_label), post_pop, "...")
                continue
            else:
                print("Creating projection {:25} "
                      "from {:10}".format(conn_label, pre_pop),
                      "to {:10}".format(post_pop),
                      "with a weight of {: 2.6f}".format(weight),
                      "uS and a delay of", delay, "ms")

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

            # TODO remove this. do better.
            if "io_to_basket" in conn_label:
                conns[conns < 0] = 2 ** 32 - 1
            else:
                assert (np.max(conns[:, 0]) < self.number_of_neurons[pre_pop]), \
                    "pre id max: {} vs. {} for {}".format(np.max(conns[:, 0]), self.number_of_neurons[pre_pop], conn_label)
                assert (np.min(conns[:, 0]).astype(int) >= 0), \
                    "pre id min: {} vs. {} for {}".format(np.min(conns[:, 0]), 0, conn_label)
                assert (np.max(conns[:, 1]) < self.number_of_neurons[post_pop]), \
                    "post id max: {} vs. {} for {}".format(np.max(conns[:, 0]), self.number_of_neurons[post_pop], conn_label)
                assert (np.min(conns[:, 1]).astype(int) >= 0), \
                    "post id min: {} vs. {} for {}".format(np.min(conns[:, 1]), 0, conn_label)
            self.connections[conn_label] = np.concatenate(
                [conns, stacked_weights, stacked_delays], axis=1)

    def build_projections(self):
        """
        Construct PyNN Projections between Populations
        """
        if self.reporting:
            # Report statistics about the populations to be built
            projection_reporting(self._raw_connectivity_info,
                                 self.number_of_neurons,
                                 self.conn_params)
        # Retrieve the Projection labels
        labels = self.connections.keys()
        # Loop over each connection in `connections`
        for conn_label in labels:
            # Retrieve saved connectivity and cast to np.ndarray
            # [:, :2] -- the columns are
            # 1. Unique G(lobal)ID (!) for pre-synaptic neuron
            # 2. Unique G(lobal)ID (!) for post-synaptic neuron
            # 3. 3D distance between somas

            pre_pop = self.conn_params[conn_label]['pre']
            post_pop = self.conn_params[conn_label]['post']
            weight = self.conn_params[conn_label]['weight']

            if post_pop in ["glomerulus", "mossy_fibers"]:
                print("Ignoring connection {:25} "
                      "terminating at".format(conn_label), post_pop, "...")
                continue

            # Adding the projection to the network
            receptor_type = "inhibitory" if weight < 0 else "excitatory"
            self.projections[conn_label] = self.sim.Projection(
                self.populations[pre_pop],  # pre-synaptic population
                self.populations[post_pop],  # post-synaptic population
                # connector includes (source, target, weight, delay)
                self.sim.FromListConnector(self.connections[conn_label]),
                receptor_type=receptor_type,  # inh or exc
                label=conn_label)  # label for connection

    def __compute_stimulus(self):
        # convert stimulation times to numpy array
        stim_times = np.asarray(self.stimulus_information['stim_times'])
        f_base = np.asarray(self.stimulus_information['f_base'])
        f_peak = np.asarray(self.stimulus_information['f_peak'])
        periodic_stimulus = self.stimulus_information['periodic_stimulus']
        percent_active = self.stimulus_information["percentage_active_fibers"] \
            if "percentage_active_fibers" in self.stimulus_information.keys() \
            else None

        no_gloms = self.number_of_neurons['glomerulus']

        print("=" * 80)
        print("Stimulation periods: ", len(stim_times),
              "lasting {} ms each".format(
                  stim_times))
        # compute number of rate changes required
        number_of_slots = int(stim_times.size)
        # compute the time at which individual rates take effect
        stim_start = np.cumsum(np.append([0], stim_times))[:number_of_slots]
        print("Thus, each period starts at ", stim_start, "ms")

        # check whether this is new style of scaffold
        # if new style --> generate activity for mf and glom
        if self.new_scaffold:
            no_mf = self.number_of_neurons['mossy_fibers']
            starts = np.ones((no_mf, number_of_slots)) * stim_start
            # compute the duration (in ms) for which individual rates are active
            durations = np.ones((no_mf, number_of_slots)) * stim_times
            # Select MFs which will fire at f_peak during stimulation
            if percent_active:
                active_mfs = np.zeros(no_mf).astype(bool)
                active_mfs[
                    np.random.choice(
                        np.arange(no_mf),
                        size=int(percent_active * no_mf),
                        replace=False)] = 1
            else:
                # all mfs will fire during stimulation
                active_mfs = np.ones(no_mf).astype(bool)
            mf_rates = np.ones((no_mf, number_of_slots)) * f_base
            mf_rates[active_mfs, 1] = f_peak
            # generate spikes for mf
            mf_spikes = create_poisson_spikes(no_mf, mf_rates,
                                              starts, durations)
            # load connectivity from mf to glom
            mf_to_glom = self.connections["mossy_to_glomerulus"]
            glom_spikes = [[] for _ in range(no_gloms)]
            if mf_to_glom is not None:
                # compute spikes for glom
                for nid, spikes in enumerate(mf_spikes):
                    # select connections corresponding to this nid as a pre neuron
                    active_connections = mf_to_glom[mf_to_glom[:, 0] == nid]
                    no_targets = active_connections.shape[0]
                    # for each row of spikes
                    for spike in spikes:
                        # for each individual spike time, post neuron and delay
                        for s, t, d in zip([spike] * no_targets,
                                           active_connections[:, 1].astype(int),
                                           active_connections[:, 3]):
                            # copy that spike with the applied delay
                            glom_spikes[t].append(s + d)

            return {'glomerulus': {'spike_times': glom_spikes},
                    'mossy_fibers': {'spike_times': mf_spikes}}
        else:  # old scaffold
            # identify ID for glomerulus population
            no_gloms = self.number_of_neurons['glomerulus']
            starts = np.ones((no_gloms, number_of_slots)) * stim_start
            # compute the duration (in ms) for which individual rates are active
            durations = np.ones((no_gloms, number_of_slots)) * stim_times
            # compute the individual rates (in Hz) for each slot
            rates = np.ones((no_gloms, number_of_slots)) * \
                    np.asarray([f_base, f_base + f_peak, f_base])
            first_glom = self.nid_offset['glomerulus']
            # extract individual gloms for the position matrix based on their
            # GID
            gloms_pos = np.empty((no_gloms, 5))
            for i in range(no_gloms):
                gloms_pos[i] = self.cell_positions[
                               self.cell_positions[:, 0] == i + first_glom, :]

            # count number of placements per id
            unique_ids = np.unique(self.cell_positions[:, 1]).astype(int)
            reverse_id_mapping = {}
            for ui in unique_ids:
                reverse_id_mapping[
                    np.count_nonzero(self.cell_positions[:, 1] == ui)] = ui
            glom_place_id = reverse_id_mapping[no_gloms]
            print("Probable placement ID for glomerulus: ", glom_place_id)
            radius = self.stimulus_information['stim_radius']
            # gloms_pos = self.cell_positions[
            #             self.cell_positions[:, 1] == glom_place_id, :]
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
            mask = np.ones(no_gloms, np.bool)
            mask[target_gloms] = 0
            # Set the firing rate of other gloms to baseline level
            rates[mask, 1] = f_base
            # Report numbers here
            print("=" * 80)
            print("Number of stimulated Gloms: ", len(target_gloms),
                  "i.e. {:6.2%} the total".format(
                      len(target_gloms) / float(no_gloms)))

            if not periodic_stimulus:
                spike_times = create_poisson_spikes(no_gloms, rates,
                                                    starts, durations)
            else:
                # SPIKE SOURCE ARRAY + PERIODIC STIMULUS
                # prepare spike times for all n_inputs neurons
                spike_times = [[] for _ in range(no_gloms)]
                for i, rate, start, duration in zip(range(no_gloms), rates, starts, durations):
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
            return {'glomerulus': {
                'spike_times': spike_times
            }}

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
            # if self.neuron_models[label] == "spikesourcearray":
            #     print("NOT enabling recordings for ", label,
            #           "(it's a Spike Source Array)")
            #     continue
            if pop:
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
            if pop:
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
            if not pop:
                continue
            if label in ["glomerulus", "mossy_fibers"]:
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
            if label in ["glomerulus", "mossy_fibers"]:
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
            except AttributeError as ae:
                print("Careful! Something happened when retrieving the "
                      "connectivity:", ae,
                      "\nRetrying using standard PyNN syntax...")
                all_connections[label] = \
                    np.array(p.get(('weight', 'delay'), format="list"))
            except TypeError as te:
                print("Connectivity is None (", te,
                      ") for connection", label)
                print("Connectivity as empty array.")
                all_connections[label] = np.array([])
        return all_connections

    def retrieve_population_names(self):
        return list(self.populations.keys())

    def retrieve_cell_params(self):
        return self.cell_params

    def retrieve_conn_params(self):
        return self.conn_params
