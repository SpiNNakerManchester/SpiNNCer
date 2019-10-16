from six import add_metaclass
from spinn_utilities.abstract_base import (
    AbstractBase, abstractproperty, abstractmethod)


@add_metaclass(AbstractBase)
class Circuit(object):

    @abstractmethod
    def get_circuit_inputs(self):
        """ Get a list of input populations belonging to the Circuit
        """

    @abstractmethod
    def get_circuit_outputs(self):
        """ Get a list of output populations belonging to the Circuit
        """
