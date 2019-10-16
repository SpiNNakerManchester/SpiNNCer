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
import numpy as np, h5py, sys
from scaffold_params import cell_type_ID
from circuit import Circuit


class Cerebellum(Circuit):
    def __init__(self, sim, connectivity):
        self.sim = sim
        self.populations = []
        self.projections = []

        self.connectivity = self.__extract_connectivity(connectivity)

        # Construct PyNN neural Populations
        self.__build_populations()

        # Construct PyNN Projections
        self.__build_projections()

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

    def __build_populations(self):
        """
        Construct PyNN Projections between Populations
        """
        pass

    def __build_projections(self):
        """
        Construct PyNN Projections between Populations
        """
        pass

    def get_circuit_inputs(self):
        pass

    def get_circuit_outputs(self):
        pass

    def get_all_populations(self):
        return self.populations



# Set up recordings



if __name__ == "__main__":
    # import sPyNNaker
    try:
        import spynnaker8 as sim
    except:
        import pyNN.spynnaker as sim
    # from spinncer_argparser import *
    import pylab as plt

    # Record simulation start time (wall clock)
    start_time = plt.datetime.datetime.now()
    connectivity_filename = 'scaffold_detailed__158.0x158.0_v3.hdf5'

    # Set up the simulation
    sim.setup(timestep=.1, min_delay=.1, max_delay=10)

    cerebellum_circuit = Cerebellum(sim, connectivity_filename)

    # Run the simulation

    # Retrieve recordings

    # Save results

    # Appropriately end the simulation
    sim.end()
