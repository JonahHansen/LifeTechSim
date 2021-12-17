# LifeTechSim
An architecture simulator for the Large Interferometer For Exoplanets (LIFE) initiative, particularly to identify optimal configurations and beam combiner designs for the nulling interferometer.

Written by Jonah Hansen, 2021 unless otherwise stated, with help from others.

See [this paper](https://www.example.com) for more details about the models in this simulation.

The engine of the simulator takes a PPop list of planets and outputs a JSON list containing the fluxes of various signal and noise components as a function of wavelengths and each observable interferometric output.

The analysis component takes these outputs and calculates the SNR, as well as plots

### Code Structure:
- analysis - folder of scripts for calculating SNR and other plots from output files
- engine - folder containing the simulator
  - engine/nullers - folder of each nulling architecture
  - engine/main_computer.py - main function for running the simulation on a single star
  - engine/planet_retrieval - convert a PPop output to a usable input for the simulator
  - engine/sim_functions - list of functions to calculate the various fluxes of each signal/noise component, among others
- jwst_backgrounds - [Git Repo from STScI](https://github.com/spacetelescope/jwst_backgrounds) for calculation of zodiacal light
- PPop - [Git Repo from kammerje](https://github.com/kammerje/P-pop) for calculating the star/planet population as an input to the simulator
- main_sim.py - script for running the simulation over all stars from PPop
- dist_dependence_sim.py - script for running the simulation on a set of simulated stars of given stellar types and distances
- other scripts - various scripts for testing

### Dependencies
- numpy
- matplotlib
- scipy
- healpy
- astropy
- cmasher
