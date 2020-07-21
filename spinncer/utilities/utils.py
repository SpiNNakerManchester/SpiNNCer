from elephant.spike_train_generation import homogeneous_poisson_process
import quantities as pq
import numpy as np
import pandas as pd
from hilbertcurve.hilbertcurve import HilbertCurve


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
    count, bin_edges = np.histogram(times, bins=bins)
    present_times_filter = count > 0
    selected_spike_times = (bin_edges[:-1])[present_times_filter]
    # Allow for multiple spikes in a timestep if that's how spike times get rounded
    rounded_spike_times = np.repeat(selected_spike_times, repeats=count[present_times_filter])
    # Check that there are the same number of spikes out as spikes in
    assert (len(rounded_spike_times) == len(times))
    return rounded_spike_times


def random_id_mapping(ids, seed=None):
    np.random.seed(seed)
    return np.random.permutation(ids)


def revert_id_mapping(ids, mapping):
    original_ids, new_ids = mapping
    argsorted_new_ids = np.argsort(new_ids)
    sorted_original_ids = original_ids[argsorted_new_ids]
    return sorted_original_ids[ids]


def apply_id_mapping(ids, mapping):
    original_ids, new_ids = mapping
    return new_ids[ids]


def revert_id_mapping_to_list(list_to_reorder, mapping):
    reordered_list = []
    original_ids, new_ids = mapping
    argsorted_new_ids = np.argsort(new_ids)
    sorted_original_ids = original_ids[argsorted_new_ids]
    for _id in sorted_original_ids:
        reordered_list.append(list_to_reorder[_id])
    return reordered_list


def apply_id_mapping_to_list(list_to_reorder, mapping):
    reordered_list = []
    for _id in mapping[1]:
        reordered_list.append(list_to_reorder[_id])
    return reordered_list


def select_cells_in_volume(positions,
                           x_min=None, x_max=None,
                           y_min=None, y_max=None,
                           z_min=None, z_max=None):
    """
    It sensible to include min or min and max in the selection?
    I'm going to just include mins here, but I need to keep an eye out for the effects of this choice.
    """
    if x_min:
        positions = positions[positions['x'] >= x_min]
    if y_min:
        positions = positions[positions['y'] >= y_min]
    if z_min:
        positions = positions[positions['z'] >= z_min]

    if x_max:
        positions = positions[positions['x'] < x_max]
    if y_max:
        positions = positions[positions['y'] < y_max]
    if z_max:
        positions = positions[positions['z'] < z_max]
    return positions


def generate_hilbert_curve(n=3, p=3):
    hc = HilbertCurve(p, n)
    npts = 2 ** (n * p)
    pts = []
    for i in range(npts):
        pts.append(hc.coordinates_from_distance(i))

    x = [pt[0] for pt in pts]
    y = [pt[1] for pt in pts]
    z = [pt[2] for pt in pts]
    return x, y, z


def hilbert_id_mapping(positions, nid_offset, no_slices=8):
    # Hilber curve
    x, y, z = generate_hilbert_curve(int(np.ceil(np.log(no_slices))))
    hilbert_df = pd.DataFrame()
    hilbert_df['x'] = x
    hilbert_df['y'] = y
    hilbert_df['z'] = z

    x_min = positions.x.min()
    x_max = positions.x.max()
    y_min = positions.y.min()
    y_max = positions.y.max()
    z_min = positions.z.min()
    z_max = positions.z.max()

    # + 0.1 because we don want to not select some cells due to value being == to max
    x_positions = np.linspace(x_min, x_max + 0.1, no_slices + 1)
    y_positions = np.linspace(y_min, y_max + 0.1, no_slices + 1)
    z_positions = np.linspace(z_min, z_max + 0.1, no_slices + 1)

    # Test the number of selected cells
    no_cells = positions.shape[0]
    rolling_total_selected_cells = 0
    no_rows = hilbert_df.shape[0]
    mapping_for_pop = np.zeros(no_cells).astype(int)
    curr_max_id = 0
    for index in range(no_rows):
        curr_row = hilbert_df.iloc[index]
        cells_in_current_voxel = select_cells_in_volume(positions,
                                                        x_min=x_positions[curr_row.x],
                                                        x_max=x_positions[curr_row.x + 1],
                                                        y_min=y_positions[curr_row.y],
                                                        y_max=y_positions[curr_row.y + 1],
                                                        z_min=z_positions[curr_row.z],
                                                        z_max=z_positions[curr_row.z + 1])
        no_cells_in_curr_voxel = cells_in_current_voxel.shape[0]
        rolling_total_selected_cells += no_cells_in_curr_voxel
        # Sort by Y
        cells_in_current_voxel = cells_in_current_voxel.sort_values(by=['y', 'z', 'x'])
        old_nids = cells_in_current_voxel['nid'].values - nid_offset
        mapping_for_pop[old_nids.astype(int)] = np.arange(curr_max_id, no_cells_in_curr_voxel + curr_max_id).astype(int)
        curr_max_id += no_cells_in_curr_voxel

    assert rolling_total_selected_cells == no_cells
    assert curr_max_id == no_cells, curr_max_id
    assert curr_max_id == np.max(mapping_for_pop) + 1, curr_max_id

    return mapping_for_pop
