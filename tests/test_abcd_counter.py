"""
This script tests whether spike counting additional provenance is correct
"""
try:
    # this might be deprecated soon
    import spynnaker8 as sim
except ImportError:
    import pyNN.spynnaker as sim
from spinncer.analysis_common import *

plt.viridis()

# ensure we use the same rc parameters for all matplotlib outputs
mlib.rcParams.update({'font.size': 24})
mlib.rcParams.update({'errorbar.capsize': 5})
mlib.rcParams.update({'figure.autolayout': True})
viridis_cmap = mlib.cm.get_cmap('viridis')
sim.setup(timestep=1, min_delay=1, max_delay=1)

# Compute 1 spike
n_neurons = 6
all_neurons = 6
spike_times = [[] for _ in range(n_neurons)]
spike_times[0] = [10, 40, 40]
spike_times[1] = [15, 45, 45]
spike_times[2] = [20, 50, 50]
spike_times[3] = [25, 55, 55]
spike_times[4] = [30, 60, 60]
spike_times[5] = [35, 65, 65]

# input pop
inp_pop = sim.Population(
    n_neurons,
    cellclass=sim.SpikeSourceArray,
    cellparams={'spike_times': spike_times},
    label="One spike pop")

# lif pop
lif_pop = sim.Population(
    all_neurons,
    cellclass=sim.IF_cond_exp,
    label="LIF pop 1 spike")

# connection
connection_list = [
    # 0 doesn't connect to anything
    [1, 0],  # 1 connects to a single post, etc
    [2, 0], [2, 1],
    [3, 0], [3, 1], [3, 2],
    [4, 0], [4, 1], [4, 2], [4, 3],
    [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5]
]
proj = sim.Projection(inp_pop, lif_pop, sim.FromListConnector(connection_list),
                      synapse_type=sim.StaticSynapse(weight=1.0, delay=1.0))
lif_pop.record(["spikes", "gsyn_exc", "v", "gsyn_inh"])
inp_pop.record(["spikes"])

sim.run(100)
val = lif_pop.get_data('v')
gsyn_exc = lif_pop.get_data('gsyn_exc')
gsyn_inh = lif_pop.get_data('gsyn_inh')
spikes = inp_pop.get_data("spikes")
out_spikes = lif_pop.get_data('spikes')
final_connectivity = \
    np.array(proj.get(('weight', 'delay'),
                      format="list")._get_data_items())

sim.end()
spikes = convert_spikes(spikes)
conn = np.array(connection_list)
inc_spike_count = np.zeros((6, 100))
per_conn_worst_spikes = []
# The following is expensive time wise
all_v = np.array(val.segments[0].filter(name='v')[0]).T
post_pop = "lif"
pre_pop = "ssa"
curr_spikes = spikes
inc_spikes_for_this_conn = np.zeros(inc_spike_count.shape)
for nid, t in curr_spikes:
    nid = int(nid)
    times = int(t)
    targets = conn[conn[:, 0] == nid][:, 1].astype(int)
    inc_spikes_for_this_conn[targets, times] += 1
inc_spike_count += inc_spikes_for_this_conn
per_conn_worst_spikes = inc_spikes_for_this_conn

nid, tstep = np.unravel_index(np.argmax(all_v, axis=None), all_v.shape)
print("{:20}-> neuron {:>8d} received {:>6d}".format(
    "LIF", int(nid), int(np.max(all_v))),
    "nA in timestep #{:8d}".format(int(tstep)))
# Also treat voltage as if it's a piggybacked value packaging
# counts of number of post-synaptic hits for spikes

all_exc_gsyn = np.array(gsyn_exc.segments[0].filter(name='gsyn_exc')[0]).T

max_v = np.max(all_v, axis=0)
max_g = np.max(all_exc_gsyn, axis=0)
abcd_view = max_v.astype(np.uint32).view(dtype=[
    #                                   ('d', np.uint16),
    #                                   ('c', np.uint16),
    ('b', np.uint16),
    ('a', np.uint16)
])
abcd_view_part_2 = max_g.astype(np.uint32).view(dtype=[
    ('d', np.uint16),
    ('c', np.uint16),
    #                                   ('b', np.uint16),
    #                                   ('a', np.uint16)
])
pd_view = pd.DataFrame(abcd_view)
pd_view_part_2 = pd.DataFrame(abcd_view_part_2)
pd_view = pd.concat([pd_view, pd_view_part_2], axis=1, sort=False)

pd_view.describe()
# this is meaningless using default tools
all_post_hits = pd_view

print("=" * 80)
print("Incoming spikes statistics")
print("-" * 80)
counts = inc_spike_count
nid, tstep = np.unravel_index(np.argmax(counts, axis=None), counts.shape)
print("{:20}-> neuron {:>8d} saw {:>6d}".format(
    "LIF", int(nid), int(np.max(counts))),
    "spikes in timestep #{:8d}".format(int(tstep)))
# Print worst case statistic for the population
maxs = np.max(counts, axis=1)
assert maxs.size == all_neurons
print("\t# spikes: mean {:8.4f}, "
      "std {:8.4f}".format(
    np.mean(maxs), np.std(maxs)
))
curr_nid_inc_spikes = per_conn_worst_spikes[nid]
print("\t\t{:10} contributed {:8d} spikes".format(
    "SSA-LIF", int(curr_nid_inc_spikes[tstep])
))

print("Plotting a b c d post-synaptic hits V2")
l = 13
transparency_lvl = .7
f = plt.figure(1, figsize=(l, l), dpi=400)
plt.close(f)

sh_exc = np.max((val.segments[0].filter(name='v')[0].magnitude.T * 2 ** 15).astype(int), axis=0)
sh_exc_2 = np.max((gsyn_exc.segments[0].filter(name='gsyn_exc')[0].magnitude.T * 2 ** 15).astype(int), axis=0)

f = plt.figure(1, figsize=(l, l), dpi=400)
a = ((sh_exc) & 0xFFFF)
b = ((sh_exc >> 16) & 0xFFFF)
c = ((sh_exc_2) & 0xFFFF)
d = ((sh_exc_2 >> 16) & 0xFFFF)

plt.plot(a + b + c + d, label='total', c="k", alpha=.3)
plt.plot(d, label='d', c="C3", alpha=transparency_lvl)
plt.plot(c, label='c', c="C2", alpha=transparency_lvl)
plt.plot(b, label='b', c="C1", alpha=transparency_lvl)
plt.plot(a, label='a', c="C0", alpha=transparency_lvl)

plt.ylabel("Count")

# plt.suptitle(use_display_name(pop))

# plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
#          stim_wanted_times.max() * time_to_bin_conversion)
# plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
plt.legend(loc="best")
plt.xlabel("Time (ms)")

plt.tight_layout()
save_figure(
    plt,
    "OR_POST_HITS_TEST",
    extensions=[".pdf", ])
plt.show()
plt.close(f)

# THIS IS BROKEN!
# print("Plotting a b c d post-synaptic hits")
# l = 13
# transparency_lvl = .7
#
# counts = all_post_hits
# f, axes = plt.subplots(2, 2, figsize=(l, l), dpi=400,
#                        sharey=True, sharex=True)
# ax_a = axes[0, 0]
# ax_b = axes[0, 1]
# ax_c = axes[1, 0]
# ax_d = axes[1, 1]
#
# ax_a.plot(counts.a, c="C0", alpha=transparency_lvl)
# ax_b.plot(counts.b, c="C1", alpha=transparency_lvl)
# ax_c.plot(counts.c, c="C2", alpha=transparency_lvl)
# ax_d.plot(counts.d, c="C3", alpha=transparency_lvl)
#
# ax_a.set_title("a")
# ax_b.set_title("b")
# ax_c.set_title("c")
# ax_d.set_title("d")
#
# ax_a.set_ylabel("# cases")
# ax_c.set_ylabel("# cases")
#
# ax_c.set_ylabel("Time (ms)")
# ax_d.set_ylabel("Time (ms)")
#
# plt.tight_layout()
# save_figure(
#     plt,
#     "POST_HITS_TEST",
#     extensions=[".pdf", ])
# plt.show()
# plt.close(f)


# Plot distribution of worst case spikes per population
print("Plotting histogram of worst spike counts")
counts = inc_spike_count
maxs = np.max(counts, axis=1)
assert maxs.size == all_neurons
f = plt.figure(1, figsize=(9, 9), dpi=400)
plt.hist(maxs, bins=20, color=viridis_cmap(0 / (1 + 1)),
         rasterized=True,
         edgecolor='k')

# plt.title(use_display_name(pop))

plt.ylabel("Count")
plt.xlabel("Max # of spikes per neuron")
plt.tight_layout()
save_figure(
    plt,
    "max_spikes_per_neuron_in_TEST",
    extensions=[".pdf", ])
plt.close(f)

print("Plotting worst_case spikes per pop")
all_ws_contribution = {'ssa_lif': per_conn_worst_spikes}
maxs = []
labels = []
for k, val in all_ws_contribution.items():
    labels.append(k)
    maxs.append(np.sum(val, axis=0))
maxs = np.asarray(maxs)

total = np.sum(maxs, axis=0)
assert total.size == maxs.shape[1]
f = plt.figure(1, figsize=(14, 9), dpi=400)
plt.plot(total, c='k', label="Total", alpha=.3)
for i, row in enumerate(maxs):
    plt.plot(row, rasterized=True, label=labels[i], alpha=.7)
# plt.title(use_display_name(pop))

# plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
#          stim_wanted_times.max() * time_to_bin_conversion)
# plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
plt.legend(loc="best")
plt.ylabel("Max # of spikes per population")
plt.xlabel("Time (ms)")
plt.tight_layout()
save_figure(
    plt,
    "sum_max_spikes_per_pop_TEST",
    extensions=[".pdf", ])
plt.close(f)

maxs = []
labels = []
for k, val in all_ws_contribution.items():
    labels.append(k)
    maxs.append(np.max(val, axis=0))
maxs = np.asarray(maxs)
total = np.sum(maxs, axis=0)
assert total.size == maxs.shape[1]
f = plt.figure(1, figsize=(14, 9), dpi=400)
plt.plot(total, c='k', label="Total", alpha=.3)
for i, row in enumerate(maxs):
    plt.plot(row, rasterized=True, label=labels[i], alpha=.7)
# plt.title(use_display_name(pop))
# stim_wanted_times = np.linspace(time_filter[1] - 50,
#                                 time_filter[2] + 50, 4).astype(int)
# plt.xlim(stim_wanted_times.min() * time_to_bin_conversion,
#          stim_wanted_times.max() * time_to_bin_conversion)
# plt.xticks(stim_wanted_times * time_to_bin_conversion, stim_wanted_times)
plt.legend(loc="best")
plt.ylabel("Max # of spikes for a neuron")
plt.xlabel("Time (ms)")
plt.tight_layout()
save_figure(
    plt,
    "worst_neuron_only_spikes_per_pop_TEST",
    extensions=[".pdf", ])
plt.close(f)

print("The end")
