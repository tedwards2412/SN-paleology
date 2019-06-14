# Paleo-detectors and neutrinos

*Is it possible find neutrinos from galactic core collapse SN with a paleo-detector?*

[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php) [![DOI](https://zenodo.org/badge/187173551.svg)](https://zenodo.org/badge/latestdoi/187173551)

Code for calculating track length spectra in paleo-detectors and exploring paleo-detector sensitivity to SN neutrinos. This code can be used to reproduce the results of [arXiv:1906.05800](http://arxiv.org/abs/1906.05800), "*Paleo-Detectors for Galactic Supernova Neutrinos*".

More information about paleo-detectors can also be found in [arXiv:1811.10549](http://arxiv.org/abs/1811.10549), [arXiv:1806.05991](http://arxiv.org/abs/1806.05991), and [arXiv:1811.06844](http://arxiv.org/abs/1811.06844).

**Authors:** Thomas D P Edwards, Bradley J Kavanagh, 

Please get in touch with any questions, comments or bug-reports.

### Overview: Paleopy

The core of the code runs on the [`paleopy`](https://github.com/tedwards2412/paleopy) package. This includes data for converting recoil energies to track lengths, along with tables of background distributions. This then allows you to calculate all the relevant track length distributions. The currently supported minerals are Nchwaningite, Sinjarite, Halite, Olivine, Gypsum and Phlogopite.

To run the notebooks it is required to install the paleopy package with:

    pip3 install git+https://github.com/tedwards2412/paleopy


### Notebooks

Details need to be added here


### Results

Tables of results can be found in  [`data/results`](data/results). The Euclideanised signals that are used to calculate many of the final results can be found in [`ES/hdf5`](ES/hdf5)

### Requirements

The code in this repo should run with Python3. Standard library requirements are in [`requirements.txt`](requirements.txt). In addition, you will need [`swordfish`](https://github.com/cweniger/swordfish) for performing the statistical analysis and [`WIMpy`](https://github.com/bradkav/WIMpy_NREFT) for calculating the DM and neutrino spectra. 
