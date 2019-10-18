import argparse

DEFAULT_FIGURE_DIR = 'figures/'
DEFAULT_RESULT_DIR = 'results/'

parser = argparse.ArgumentParser(
    description='Run a cerebellar simulation on SpiNNaker.',
    formatter_class=argparse.RawTextHelpFormatter)

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

parser.add_argument('--figures_dir', type=str,
                    help='directory into which save figures',
                    default=DEFAULT_FIGURE_DIR)

parser.add_argument('--result_dir', type=str,
                    help='directory into which to save simulation results',
                    default=DEFAULT_RESULT_DIR)

args = parser.parse_args()