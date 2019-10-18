# SpiNNCer
## A SpiNNaker-based Cerebellum simulation

This project implements a model of a Cerebellum in PyNN on SpiNNaker.

The full description of the model is available at: https://www.frontiersin.org/articles/10.3389/fninf.2019.00037


The architecture of the scaffold model:
![Cerebellum architecture](https://www.frontiersin.org/files/Articles/444802/fninf-13-00037-HTML-r1/image_m/fninf-13-00037-g001.jpg)

Neuron parameters:
![Neuron parameter table](https://www.frontiersin.org/files/Articles/444802/fninf-13-00037-HTML-r1/image_m/fninf-13-00037-t002.jpg)

## Installation

In a (virtualenv) console type the usual:

`python setup.py develop`, optionally with `--no-deps`

or

`pip install .`


## Use

[spinncer/cerebellum_on_spinnaker.py](spinncer/cerebellum_on_spinnaker.py) 
is the main script of the repository. It implements the `Cerebellum` class.
Instantiating this class (i.e. `cerebellum = Cerebellum(...)`) has to be done
**AFTER** calling the PyNN `sim.setup(...)` construct. It will  automatically
add the appropriate elements to the PyNN network 
(i.e. `Populations` and `Projections`).

Additional information is available by typing:
`python cerebellum_on_spinnaker.py -h`



