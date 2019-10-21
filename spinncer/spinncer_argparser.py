import argparse

DEFAULT_FIGURE_DIR = 'figures/'
DEFAULT_RESULT_DIR = 'results/'
DEFAULT_SIMTIME = 500  # ms
DEFAULT_TIMESTEP = 0.1  # ms

parser = argparse.ArgumentParser(
    description='Run a cerebellar simulation on SpiNNaker.',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--simtime', type=float,
                    help="simulation time (in ms) "
                         "-- default {}ms".format(DEFAULT_SIMTIME),
                    default=DEFAULT_SIMTIME)


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