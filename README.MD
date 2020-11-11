# SpiNNCer
## Neuromorphic cerebellum implementation on SpiNNaker

This work presents the first simulation of a large-scale, bio-inspired cerebellum model performed on neuromorphic hardware. A network containing 97k neurons and 4.2M synapses is simulated on the SpiNNaker neuromorphic system. The same model is also simulated using NEST, a popular spiking neural network simulator using generic computational resources and double precision floating point arithmetic, as a basis for validation. Individual cell and network-level spiking activity is validated against NEST and SpiNNaker is shown to produce results in agreement with NEST. The goal after validation is understanding how to accelerate the simulation speed of the network on the SpiNNaker system. Through detailed communication profiling, peak network activity is identified as one of the main challenges for simulation speed-up.

Supporting data for this project is available at: https://kg.ebrains.eu/search/instances/Dataset/069e4718-2b4a-4056-9c8f-e646a841b2dd

The full description of the model is available at: https://www.frontiersin.org/articles/10.3389/fninf.2019.00037

## Installation


The software package can be installed in a local virtual environment. After unzipping SpiNNCer-master.zip, the following can be run inside the unarchived folder to install the package:


`python setup.py develop`

or

`pip install .`

This installation assumes the existence of either sPyNNaker or NEST, for brevity. For installation instructions for those, please follow the official directions provided by each package.


## Use


After installation, simulations can be run in the spinncer directory e.g. the following will generate the reproducible results:
```
python cerebellum_experiment.py --simulator nest --input scaffold_full_dcn_400.0x400.0_v3.hdf5 -o nest_400x400_pss_3 -s 400x400_stimulus_3.npz --id_remap grid
```
Where `scaffold_full_dcn_400.0x400.0_v3.hdf5` is the cerebellar model to be instantiated in PyNN,   `-o nest_400x400_pss_3` the name of the output file, 
`-s 400x400_stimulus_3.npz` the location of a stimulus file if desired.

For SpiNNaker, and additional argument may be required depending on whether weight and injected current normalisation is enabled:

```
python cerebellum_experiment.py --simulator spinnaker --input scaffold_full_dcn_400.0x400.0_v3.hdf5 -o spinnaker_400x400_pss_3 -s 400x400_stimulus_3.npz --id_remap grid --r_mem
```

**For further information run:**
```
python cerebellum_experiment.py --help
python cerebellum_analysis.py --help
```
Analysing a large-scale experiment
Analysing the results of individual experiments can be performed by running:
```
python cerebellum_analysis.py --consider_delays worst_case_spike   --highlight_stim -i results/spinn_400x400_pss_3.npz 
```
This will produce both textual output (which can be piped to a file) and figures (these will be placed in the figures/ directory). 
`--consider_delays` informs the analysis to adjust time windows for individual populations in the model to account for the delay in spike propagation defined by the prescribed synaptic delays. `--worst-case-spike` informs the analysis to report information about peak activity in the network. `--highlight_stim` informs the analysis to produce figures in which the stimulation period is highlighted by shading the appropriate time period.
Comparing two large-scale experiments
Comparing two simulations can be performed by running:
```
python cerebellum_analysis.py --consider_delays --compare results/spinn_400x400_pss_3.npz  results/nest_400x400_pss_3.npz
```
This will produce both textual output (which can be redirected to a file) and figures (these will be placed in the figures/ directory). 
`--consider_delays` informs the analysis to adjust time windows for individual populations in the model to account for the delay in spike propagation defined by the prescribed synaptic delays.
Analysing multiple runs of the large-scale experiment
Analysing multiple experiments is used mainly to profile the peak number of spikes received per core on SpiNNaker. This analysis can be run by pointing the provenance_analysis.py to the directory (or directories in this case) containing the multiple simulations to be analysed. The directories are analysed in order. 
```
python provenance_analysis.py -i activity_sweep_f_peak_POISSON_@1000x_vanilla_scaling_200_grid activity_sweep_f_peak_PERIODIC_@1000x_vanilla_scaling_200_grid        --group_on f_peak
python provenance_analysis.py -i activity_sweep_stim_radius_POISSON_@1000x_vanilla_scaling_200_grid activity_sweep_stim_radius_PERIODIC_@1000x_vanilla_scaling_200_grid        --group_on stim_radius
```
Here, `--group_on` informs the analysis of which parameter is being varied in the set of experiments to be analysed.

[spinncer/cerebellum.py](spinncer/cerebellum.py) implements the `Cerebellum` class.
Instantiating this class (e.g. `cerebellum = Cerebellum(...)`) has to be done
**AFTER** calling the PyNN `sim.setup(...)` construct. It will  automatically
add the appropriate elements to the PyNN network 
(i.e. `Populations` and `Projections`).



