import argparse

DEFAULT_FIGURE_DIR = 'figures/'

analysis_parser = argparse.ArgumentParser(
    description='Analyse a cerebellar simulation run on SpiNNaker or  NEST.',
    formatter_class=argparse.RawTextHelpFormatter)

analysis_parser.add_argument('-i', '--input', type=str, nargs="*",
                             help="name(s) of the npz archive storing "
                                  "the results from running the simulations",
                             dest='input')

analysis_parser.add_argument('--figures_dir', type=str,
                             help='directory into which to save figures',
                             default=DEFAULT_FIGURE_DIR)
analysis_args = analysis_parser.parse_args()
