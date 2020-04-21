# InhPlastVisualization
Simulations for Inhibitory Plasticity visualization tools. 

This code corresponds to the paper: 
Simulation, visualization and analysis tools for pattern recognition assessment with spiking neuronal networks
by Sergio E. Galindo, Pablo Toharia, Óscar D. Robles, Eduardo Ros, Luis Pastor and Jesús A. Garrido

Please, if you work in research and publish a paper, do not hesitate to cite this publication as follows:
Galindo, S. E., Toharia, P., Robles, O. D., Ros, E., Pastor, L., & Garrido, J. A. (2020). Simulation, visualization and analysis tools for pattern recognition assessment with spiking neuronal networks. Neurocomputing.

You can find the original article in this link: https://doi.org/10.1016/j.neucom.2020.02.114, and you 
may also find a post-print version of the article here: https://arxiv.org/abs/2003.06343

This repository includes the source code implemented for the simulation of the spiking neural network with
excitatory and inhibitory plasticity, the optimization with an evolutionary algorithm and the analysis of the data.

Please, note that the source code of the main visualization tool, ViSimpl can be found in its own repository: https://github.com/gmrvvis/visimpl.

## Dependences
Most of the code in this repository has been implemented in Python and the following packages can be required:
- Nest (2.14 has been used to run the simulations in this repository).
- Numpy
- SciPy
- Matplotlib
- Mpi4py

## Folder structure
The following directories can be found in this repository:
- src: It includes the python source code to simulate the neuronal networks and the evolutionary algorithms.
- launch/GoCFanIn: It includes all the scripts required to launch the Evolutionary Algorithm optimization. Please,
note that the evolutionary algorithms were launched in parallel in a cluster with ~500 nodes available running for several days,
so it might take a long time.
- config/GoCFanIn: It includes the configuration files for the Evolutionary Algorithms and the neuronal networks.
- analysis/GoCFanIn: It includes Jupyter notebooks for analyzing the results from the EA and running the best-found configurations.
- results/GoCFanIn: It includes the resulting data of the Evolutionary Algorithms as well as the log files of these simulations.

If you have any additional question/comment/suggestion, please do not hesitate to contact me.

Jesús Garrido
University of Granada
