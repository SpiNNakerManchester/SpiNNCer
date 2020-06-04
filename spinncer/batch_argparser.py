import argparse

# Current defaults for [App: Motion detection]
# as of 08.10.2018

# Default values
DEFAULT_NO_CPUS = 32
DEFAULT_MAX_CONCURRENT_PROCESSES = 24
DEFAULT_SUFFIX = None

# Argument parser
parser = argparse.ArgumentParser(
    description='Batch runner adapted from '
                'https://github.com/pabogdan/neurogenesis/blob/master/synaptogenesis/batch_argparser.py'
                'and '
                'https://github.com/pabogdan/neurogenesis/blob/master/synaptogenesis/batch_runner.py')

parser.add_argument('path', help='path of input .npz archive defining '
                                 'connectivity', nargs='*')

parser.add_argument('-o', '--output', type=str,
                    help="name of the numpy archive storing simulation results",
                    dest='filename')

parser.add_argument('--suffix', type=str,
                    help="add a recognisable suffix to all the file "
                         "generated in this batch "
                         "-- [default {}]".format(DEFAULT_SUFFIX),
                    dest='suffix')

parser.add_argument('--no_cpus', type=int,
                    default=DEFAULT_NO_CPUS, dest='no_cpus',
                    help='total number of available CPUs'
                         ' -- [default {}]'.format(DEFAULT_NO_CPUS))

parser.add_argument('--max_processes', type=int,
                    default=DEFAULT_MAX_CONCURRENT_PROCESSES,
                    dest='max_processes',
                    help='total number of concurrent processes'
                         ' -- [default {}]'.format(DEFAULT_MAX_CONCURRENT_PROCESSES))

args = parser.parse_args()
