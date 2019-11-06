import argparse

DEFAULT_FIGURE_DIR = 'figures/'
DEFAULT_RESULT_DIR = 'results/'
DEFAULT_SIMTIME = 1000  # ms, see page 5 of 10.3389/fninf.2019.00037
DEFAULT_TIMESTEP = 0.1  # ms

# Stimulus: 300 ms pre-stimulus, 50 ms stimulus, and 650 ms post-stimulus.
DEFAULT_STIMULATION_DURATIONS = [300, 50, 650]
DEFAULT_BACKGROUND_RATE = 1.  # Hz
DEFAULT_BURST_RATE = 150.  # Hz

# Radius of volume of excited gloms
STIMULATED_GLOMERULUS_RADIUS = 140  # micrometers

parser = argparse.ArgumentParser(
    description='Run a cerebellar simulation on SpiNNaker.',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--simtime', type=float,
                    help="simulation time (in ms) "
                         "-- default {}ms".format(DEFAULT_SIMTIME),
                    default=DEFAULT_SIMTIME)

parser.add_argument('--periodic_stimulus', action="store_true",
                    help='if this flag is present the stimulus '
                         'is a spike source array with periodic stimulus')

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

parser.add_argument('-o', '--output', type=str,
                    help="name of the numpy archive "
                         "storing simulation results",
                    dest='filename')

parser.add_argument('--no_reports', action="store_false",
                    help='if this flag is present suppress the '
                         'printing of textual reports')

parser.add_argument('--skip_projections', action="store_true",
                    help='if this flag is present suppress the '
                         'creation of projections in the network')

parser.add_argument('--figures_dir', type=str,
                    help='directory into which save figures',
                    default=DEFAULT_FIGURE_DIR)

parser.add_argument('--result_dir', type=str,
                    help='directory into which to save simulation results',
                    default=DEFAULT_RESULT_DIR)

args = parser.parse_args()
