import argparse

DEFAULT_FIGURE_DIR = 'figures/'
DEFAULT_RESULT_DIR = 'results/'

# Default dataset file
DEFAULT_DATASET = 'scaffold_detailed__158.0x158.0_v3.hdf5'


parser = argparse.ArgumentParser(
    description='Analyse a cerebellar simulation run on SpiNNaker or  NEST.',
    formatter_class=argparse.RawTextHelpFormatter)


parser.add_argument('-i', '--input', type=str, nargs="*",
                    help="name of the npz archive storing "
                         "the results from running the simulations",
                    dest='input')

parser.add_argument('--figures_dir', type=str,
                    help='directory into which to save figures',
                    default=DEFAULT_FIGURE_DIR)
args = parser.parse_args()
