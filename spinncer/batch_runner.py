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

# Current defaults for [App: Motion detection]
# as of 26.09.2018

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
PHASES = [POISSON_PHASE, PERIODIC_PHASE]
PHASES_NAMES = ["poisson", "periodic"]
PHASES_ARGS = [None, "--periodic_stimulus"]

concurrently_active_processes = 0

rb_left_shifts = np.arange(16)

# Compute total number of runs
total_runs = rb_left_shifts.size * len(PHASES)

parameters_of_interest = {
    'rb_left_shifts': rb_left_shifts,
    'phase': PHASES,
}

dataset = "scaffold_full_dcn_400.0x400.0_v3.hdf5"

log_calls = []

# making a directory for this experiment
dir_name = "ls_sweep_@{}".format(suffix)
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

for phase in PHASES:
    for ls in rb_left_shifts:
        curr_params = (ls, phase)
        filename = "spinn_400x400" \
                   "_ls_{}" \
                   "_{}" \
                   "_@{}".format(ls,
                                 PHASES_NAMES[phase],
                                 suffix)
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

        if not prev_run:
            os.mknod("results.etxt")
            os.mknod("results.otxt")

        concurrently_active_processes += 1
        null = open(os.devnull, 'w')
        print("Run ", concurrently_active_processes, "...")

        call = [sys.executable,
                '../../cerebellum_experiment.py',
                '--input', dataset,
                '-o', filename,
                '--rb_left_shifts'
                ]
        for _ in range(2):
            call.append(str(ls))

        if PHASES_ARGS[phase] is not None:
            call.append(PHASES_ARGS[phase])
        print("CALL", call)
        log_calls.append((call, filename, curr_params))
        if (concurrently_active_processes % MAX_CONCURRENT_PROCESSES == 0
                or concurrently_active_processes == total_runs):
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
        # TODO block if re-running simulations and not yet done (how would I know?)
print("All done!")

end_time = plt.datetime.datetime.now()
total_time = end_time - currrent_time
np.savez_compressed("batch_{}".format(suffix),
                    parameters_of_interest=parameters_of_interest,
                    total_time=total_time,
                    log_calls=log_calls)
