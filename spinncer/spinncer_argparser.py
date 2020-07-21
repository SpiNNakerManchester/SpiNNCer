import argparse

DEFAULT_FIGURE_DIR = 'figures/'
DEFAULT_RESULT_DIR = 'results/'
DEFAULT_SIMTIME = 1000  # ms, see page 5 of 10.3389/fninf.2019.00037
DEFAULT_TIMESTEP = 0.1  # ms

# Stimulus: 300 ms pre-stimulus, 50 ms stimulus, and 650 ms post-stimulus.
DEFAULT_STIMULATION_DURATIONS = [300, 50, 650]
DEFAULT_BACKGROUND_RATE = 1.  # Hz
DEFAULT_BURST_RATE = 150.  # Hz

# Simulator
DEFAULT_SIMULATOR = "spinnaker"

# Neuron model
DEFAULT_NEURON_MODEL = "IF_cond_exp"

# Scale weights from the default
DEFAULT_WEIGHT_SCALING = 1.

# Default dataset file
DEFAULT_DATASET = 'scaffold_detailed__158.0x158.0_v3.hdf5'
DEFAULT_PARAMETER_JSON = None

# Radius of volume of excited gloms -- Paper misreports this number,
# it's actually 70 micrometers for the 400x400 um2 model
STIMULATED_GLOMERULUS_RADIUS = 140  # micrometers

DEFAULT_PERCENTAGE_ACTIVE_FIBERS = .8

DEFAULT_TIMESCALE = None

DEFAULT_RB_LEFT_SHIFTS = None

DEFAULT_NO_LOOPS_GRC = 3

DEFAULT_SINGLE_CELL_TEST = 'granule'

parser = argparse.ArgumentParser(
    description='Run a cerebellar simulation written in PyNN '
                'using sPyNNaker on SpiNNaker or using NEST on (H)PC.',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--simtime', type=float,
                    help="simulation time (in ms) "
                         "-- default {} ms".format(DEFAULT_SIMTIME),
                    default=DEFAULT_SIMTIME)

parser.add_argument('-n', '--neuron_model', type=str,
                    help="point neuron model to use in the simulation "
                         "-- default {}".format(DEFAULT_NEURON_MODEL),
                    default=DEFAULT_NEURON_MODEL)

parser.add_argument('--weight_scaling', type=float,
                    help="scale the weights before passing them to sPyNNaker "
                         "-- default {}".format(DEFAULT_WEIGHT_SCALING),
                    default=DEFAULT_WEIGHT_SCALING)

parser.add_argument('--percentage_active_fibers', type=float,
                    help="percentage of active fibers during stimulation "
                         "-- default {}".format(DEFAULT_PERCENTAGE_ACTIVE_FIBERS),
                    default=DEFAULT_PERCENTAGE_ACTIVE_FIBERS)

parser.add_argument('--simulator', type=str,
                    help="which simulator to use "
                         "-- default {} ms".format(DEFAULT_SIMULATOR),
                    default=DEFAULT_SIMULATOR)

parser.add_argument('--periodic_stimulus', action="store_true",
                    help='if this flag is present the stimulus '
                         'is a spike source array with periodic stimulus')

parser.add_argument('-g', '--generate_conversion_constants', action="store_true",
                    help='if this flag is present the simulation is not run. '
                         'Instead, we generate a "conversion_constants.npz" '
                         'file to reduce chance of incorrectly translated '
                         'spikes from group of .npy files.')

parser.add_argument('--stim_radius', type=float,
                    help="stimulation volume radius "
                         "-- default {} um".format(STIMULATED_GLOMERULUS_RADIUS),
                    default=STIMULATED_GLOMERULUS_RADIUS)

parser.add_argument('--f_base', type=float,
                    help="background firing rate of the stimulus "
                         "-- default {}Hz".format(DEFAULT_BACKGROUND_RATE),
                    default=DEFAULT_BACKGROUND_RATE)

parser.add_argument('--f_peak', type=float,
                    help="burst firing rate of the stimulus "
                         "-- default {}Hz".format(DEFAULT_BURST_RATE),
                    default=DEFAULT_BURST_RATE)

parser.add_argument('--stim_times', type=float, nargs="+",
                    help="stimulations times (in ms) "
                         "-- default {}ms".format(DEFAULT_STIMULATION_DURATIONS),
                    default=DEFAULT_STIMULATION_DURATIONS)

parser.add_argument('--timestep', type=float,
                    help="simulation timestep (in ms) "
                         "-- default {}ms".format(DEFAULT_TIMESTEP),
                    default=DEFAULT_TIMESTEP)

parser.add_argument('-i', '--input', type=str,
                    help="name of the dataset storing "
                         "initial connectivity for the simulation",
                    dest='dataset')

parser.add_argument('--param_json', type=str,
                    help="name of the json storing "
                         "cell and connection parameters"
                         "-- default {}".format(DEFAULT_PARAMETER_JSON))

parser.add_argument('-o', '--output', type=str,
                    help="name of the numpy archive (.npz) "
                         "storing simulation results",
                    dest='filename')

parser.add_argument('--no_reports', action="store_false",
                    help='if this flag is present suppress the '
                         'printing of textual reports')

parser.add_argument('--skip_projections', action="store_true",
                    help='if this flag is present suppress the '
                         'creation of projections in the network')

parser.add_argument('--disable_around', action="store_false",
                    help='if this flag is present disable the rounding of '
                         'input spike times to nearest time step')

parser.add_argument('--figures_dir', type=str,
                    help='directory into which to save figures',
                    default=DEFAULT_FIGURE_DIR)

parser.add_argument('--result_dir', type=str,
                    help='directory into which to save simulation results',
                    default=DEFAULT_RESULT_DIR)

parser.add_argument('--rb_left_shifts', type=float, nargs="+",
                    help="left shifts for ring buffer values "
                         "-- default {}ms".format(DEFAULT_RB_LEFT_SHIFTS),
                    default=DEFAULT_RB_LEFT_SHIFTS)

parser.add_argument('--timescale', type=int,
                    help='set the slowdown factor of the simulation',
                    default=DEFAULT_TIMESCALE)

parser.add_argument('--loops_grc', type=int,
                    help='number of times around the neuron update loop for grc',
                    default=DEFAULT_NO_LOOPS_GRC)

parser.add_argument('-s', '--stimulus', type=str,
                    help="name of the file storing "
                         "input spikes",
                    dest='stimulus_from_file')

parser.add_argument('--single_spike_exp', action="store_true",
                    help='if this flag is present the input consists of a '
                         'single spike at time = 10 ms')

parser.add_argument('--disable_i_offset', action="store_true",
                    help='if this flag is present the input consists of a '
                         'single spike at time = 10 ms')

parser.add_argument('--r_mem', action="store_true",
                    help='pre-multiply membrane resistance into weight and '
                         'i_offset')

parser.add_argument('--test_max_weight', action="store_true",
                    help='[testing_all_connections.py] use the maximum weight '
                         'of expected spikes for the test (send 1 spike, scale up weight)')

parser.add_argument('--test_max_spikes', action="store_true",
                    help='[testing_all_connections.py] use the maximum number '
                         'of expected spikes for the test')

parser.add_argument('--population', type=str,
                    help="name of the population to test single cell from",
                    default=DEFAULT_SINGLE_CELL_TEST)

parser.add_argument('--id_remap', type=str,
                    help="how to remap neuron IDs",
                    default=None)

parser.add_argument('--spike_seed', type=int,
                    help='seed used for spike generation',
                    default=None)

parser.add_argument('--id_seed', type=int,
                    help='seed used for id remapping',
                    default=None)

args = parser.parse_args()
