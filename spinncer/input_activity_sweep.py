"""
Batch runner adapted from '
'https://github.com/pabogdan/neurogenesis/blob/master/synaptogenesis/batch_argparser.py'
'and '
'https://github.com/pabogdan/neurogenesis/blob/master/synaptogenesis/batch_runner.py
"""

import subprocess
import os
import numpy as np
import sys
import hashlib
import pylab as plt
from spinncer.batch_argparser import *
import shutil
from collections import namedtuple
import spinncer as x

base_path = os.path.dirname(x.__file__)

currrent_time = plt.datetime.datetime.now()
string_time = currrent_time.strftime("%H%M%S_%d%m%Y")

if args.suffix:
    suffix = args.suffix
else:
    suffix = hashlib.md5(string_time.encode('utf-8')).hexdigest()

# Some constants
NO_CPUS = args.no_cpus
MAX_CONCURRENT_PROCESSES = args.max_processes

POISSON_PHASE = 0
PERIODIC_PHASE = 1
PHASES_NAMES = ["poisson", "periodic"]
PHASES_ARGS = [None, "--periodic_stimulus"]
# Both phases
# PHASES = [POISSON_PHASE, PERIODIC_PHASE]


# PHASES = [POISSON_PHASE] # Only Poisson phase
PHASES = [PERIODIC_PHASE]  # Only PERIODIC phase

concurrently_active_processes = 0
Result = namedtuple('Result', 'call filename parameters')

# f_peaks = np.arange(30, 200, 20)  # Hz
f_peaks = np.asarray([150])  # Hz
# radii = np.asarray([140])  # um
radii = np.arange(20, 190, 20)  # um
# TODO Fix ring buffer left shift value from previous experiment
RB_LEFT_SHIFT = None

# Compute total number of runs
total_runs = f_peaks.size * len(PHASES) * radii.size

parameters_of_interest = {
    'f_peak': f_peaks,
    'stim_radius': radii,
    'phase': PHASES,
}

dataset = "scaffold_full_dcn_400.0x400.0_v3.hdf5"

log_calls = []

# making a directory for this experiment
# dir_name = "activity_sweep_stim_radius_3_loops_10x_PERIODIC_@{}".format(suffix)
# dir_name = "activity_sweep_f_peak_POISSON_@{}".format(suffix)
dir_name = "activity_sweep_stim_radius_PERIODIC_@{}".format(suffix)

print("=" * 80)
print("TOTAL RUNS", total_runs)
if not os.path.isdir(dir_name):
    print("MKDIR", dir_name)
    os.mkdir(dir_name)
else:
    print("FOLDER ALREADY EXISTS. RE-RUNNING INCOMPLETE JOBS.")
print("CHDIR", dir_name)
os.chdir(dir_name)
print("GETCWD", os.getcwd())
print("-" * 80)

params = {}
no_runs = np.arange(10)

for phase in PHASES:
    for n_run in no_runs:
        for f_peak in f_peaks:
            for stim_radius in radii:
                curr_params = {'stim_radius': stim_radius,
                               'f_peak': f_peak,
                               'n_run': n_run,
                               'phase': phase}
                filename = "spinn_400x400" \
                           "_f_peak_{}" \
                           "_run_{}" \
                           "_stim_radius_{}" \
                           "_{}" \
                           "_@{}".format(f_peak,
                                         n_run,
                                         stim_radius,
                                         PHASES_NAMES[phase],
                                         suffix)

                params[filename] = curr_params

                concurrently_active_processes += 1
                null = open(os.devnull, 'w')
                print("Run ", concurrently_active_processes, "...")

                call = [sys.executable,
                        os.path.join(base_path, 'cerebellum_experiment.py'),
                        '--input', dataset,
                        '-o', filename,
                        '--f_peak', str(f_peak),
                        '--stim_radius', str(stim_radius),
                        '--id_remap', 'grid'
                        ]

                if PHASES_ARGS[phase] is not None:
                    call.append(PHASES_ARGS[phase])

                print("CALL", call)
                log_calls.append(Result(call, filename, curr_params))

                # making a directory for this individual experiment
                prev_run = True
                if os.path.isdir(filename) and os.path.isfile(
                        os.path.join(filename, "structured_provenance.csv")):
                    print("Skipping", filename)
                    continue
                elif not os.path.isdir(filename):
                    os.mkdir(filename)
                    prev_run = False
                os.chdir(filename)
                print("GETCWD", os.getcwd())
                shutil.copyfile("../../spynnaker.cfg", "spynnaker.cfg")
                if concurrently_active_processes % MAX_CONCURRENT_PROCESSES == 0 \
                        or concurrently_active_processes == total_runs:
                    # Blocking
                    with open("results.out", "wb") as out, open("results.err", "wb") as err:
                        subprocess.call(call, stdout=out, stderr=err)
                    print("{} sims done".format(concurrently_active_processes))
                else:
                    # Non-blocking
                    with open("results.out", "wb") as out, open("results.err", "wb") as err:
                        subprocess.Popen(call, stdout=out, stderr=err)
                os.chdir("..")
                print("=" * 80)
print("All done!")

end_time = plt.datetime.datetime.now()
total_time = end_time - currrent_time
np.savez_compressed("batch_{}".format(suffix),
                    parameters_of_interest=parameters_of_interest,
                    total_time=total_time,
                    log_calls=log_calls,
                    params=params)
