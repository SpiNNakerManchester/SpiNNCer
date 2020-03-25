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

analysis_parser.add_argument('--worst_case_spikes', action="store_true",
                             help='if this flag is present the expensive '
                                  'process of counting number of afferent '
                                  'spikes is performed.')

analysis_parser.add_argument('--consider_delays', action="store_true",
                             help='if this flag is present delays in the '
                                  'propagation of excitation are taken into '
                                  'account when report firing rates.')

analysis_parser.add_argument('--dark_background', action="store_true",
                             help='if this flag is set plots will have a '
                                  'black background instead of white.')

analysis_parser.add_argument('--highlight_stim', action="store_true",
                             help='if this flag is set plots will highlight'
                                  'the stimulus time by using a different '
                                  'background colour.')

analysis_parser.add_argument('--group_on', type=str, nargs="+",
                             help="For certain kinds of analysis (e.g. "
                                  "provenance) group results based on a "
                                  "parameter of interest",
                             default=None)

analysis_parser.add_argument('--group_on_name', type=str, nargs="+",
                             help="For certain kinds of analysis (e.g. "
                                  "provenance) group results based on string "
                                  "pattern matching",
                             default=None)

analysis_args = analysis_parser.parse_args()
