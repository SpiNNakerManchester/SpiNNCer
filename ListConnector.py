from __future__ import division
from pyNN import errors
from pyNN.standardmodels import StandardSynapseType
import numpy
from pyNN.connectors import *

try:
    from itertools import izip
except ImportError:  # python3.x
    izip = zip
try:
    basestring
except NameError:
    basestring = str

import logging
from copy import deepcopy


try:
    import csa
    haveCSA = True
except ImportError:
    haveCSA = False

logger = logging.getLogger("PyNN")

class ListConnector(Connector):
    """
    Make connections according to a list.
    Arguments:
        `conn_list`:
            a list of tuples, one tuple for each connection. Each tuple should contain:
            `(pre_idx, post_idx, p1, p2, ..., pn)` where `pre_idx` is the index
            (i.e. order in the Population, not the ID) of the presynaptic
            neuron, `post_idx` is the index of the postsynaptic neuron, and
            p1, p2, etc. are the synaptic parameters (e.g. weight, delay,
            plasticity parameters).
        `column_names`:
            the names of the parameters p1, p2, etc. If not provided, it is
            assumed the parameters are 'weight', 'delay' (for backwards
            compatibility). This should be specified using a tuple.
        `safe`:
            if True, check that weights and delays have valid values. If False,
            this check is skipped.
        `callback`:
            if True, display a progress bar on the terminal.
    """
    parameter_names = ('conn_list',)

    def __init__(self, conn_list, column_names=None, safe=True, callback=None):
        """
        Create a new connector.
        """
        Connector.__init__(self, safe=safe, callback=callback)
        self.conn_list = numpy.array(conn_list)
        if len(conn_list) > 0:
            n_columns = self.conn_list.shape[1]
            if column_names is None:
                if n_columns == 2:
                    self.column_names = ()
                elif n_columns == 4:
                    self.column_names = ('weight', 'delay')
                else:
                    raise TypeError("Argument 'column_names' is required.")
            else:
                self.column_names = column_names
                if n_columns != len(self.column_names) + 2:
                    raise ValueError("connection list has %d parameter columns, but %d column names provided." % (
                                    n_columns - 2, len(self.column_names)))
        else:
            self.column_names = ()


    def connect(self, projection):
        """Connect-up a Projection."""
        logger.debug("conn_list (original) = \n%s", self.conn_list)
        synapse_parameter_names = projection.synapse_type.get_parameter_names()
        for name in self.column_names:
            if name not in synapse_parameter_names:
                raise ValueError("%s is not a valid parameter for %s" % (
                                 name, projection.synapse_type.__class__.__name__))
        if self.conn_list.size == 0:
            return
        #if numpy.any(self.conn_list[:, 0] >= projection.pre.size):
         #   raise errors.ConnectionError("source index out of range")
        # need to do some profiling, to figure out the best way to do this:
        #  - order of sorting/filtering by local
        #  - use numpy.unique, or just do in1d(self.conn_list)?
        idx = numpy.argsort(self.conn_list[:, 1])
        targets = numpy.unique(self.conn_list[:, 1]).astype(numpy.int)
        local = numpy.in1d(targets,
                           numpy.arange(projection.post.size)[projection.post._mask_local],
                           assume_unique=True)
        local_targets = targets[local]
        self.conn_list = self.conn_list[idx]
        left = numpy.searchsorted(self.conn_list[:, 1], local_targets, 'left')
        right = numpy.searchsorted(self.conn_list[:, 1], local_targets, 'right')
        logger.debug("idx = %s", idx)
        logger.debug("targets = %s", targets)
        logger.debug("local_targets = %s", local_targets)
        logger.debug("conn_list (sorted by target) = \n%s", self.conn_list)
        logger.debug("left = %s", left)
        logger.debug("right = %s", right)
        sources = []
        for tgt, l, r in zip(local_targets, left, right):
            #for i in self.conn_list[l:r, 0]:
                #print(i[0])
             #   sources.append(int(i[0]))
            sources = self.conn_list[l:r, 0].astype(numpy.int)

            connection_parameters = deepcopy(projection.synapse_type.parameter_space)
            connection_parameters.shape = (r - l,)
            for col, name in enumerate(self.column_names, 2):
                #print(col)
                #print(name)
                #print(self.conn_list[l:r, col])
                connection_parameters.update(**{name: self.conn_list[l:r, col]})
            if isinstance(projection.synapse_type, StandardSynapseType):
                connection_parameters = projection.synapse_type.translate(
                                            connection_parameters)
            connection_parameters.evaluate()
            projection._convergent_connect(sources, tgt, **connection_parameters)
