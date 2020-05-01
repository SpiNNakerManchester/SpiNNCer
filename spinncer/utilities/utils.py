from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq
import numpy as np


# Originally from, with modification
# https://codereview.stackexchange.com/questions/21033/flatten-dictionary-in-python-functional-style
def flatten_dict(d):
    def expand(key, value):
        if isinstance(value, dict):
            return [(k, v) for k, v in flatten_dict(value).items()]
        else:
            return [(key, value)]

    items = [item for k, v in d.items() for item in expand(k, v)]

    return dict(items)


def create_poisson_spikes(n_inputs, rates, starts, durations):
    spike_times = [[] for _ in range(n_inputs)]
    for i, rate, start, duration in zip(range(n_inputs), rates, starts, durations):
        curr_spikes = []
        for r, s, d in zip(rate, start, duration):
            curr_spikes.append(homogeneous_poisson_process(
                rate=r * pq.Hz,
                t_start=s * pq.ms,
                t_stop=(s + d) * pq.ms,
                as_array=True))
        spike_times[i] = np.concatenate(curr_spikes)
    return spike_times


# Originally wanted to use something like this:
# https://stackoverflow.com/questions/58065055/floor-and-ceil-with-number-of-decimals
# but the solution is wrong: e.g.
# x = 0.5 * pq.ms
# my_floor(x, 1)
# result: array(0.4) * ms
# def floor_spike_time(a, precision=0, dtype=1 * pq.ms):
#     return np.round(a - 0.5 * 10**(-precision) * dtype, precision)
def floor_spike_time(times, dt=0.1, t_start=0, t_stop=1000.0):
    bins = np.arange(t_start, t_stop + dt, dt)
    count, bin_edges = np.histogram(times , bins=bins)
    return (bin_edges[:-1])[count > 0]

