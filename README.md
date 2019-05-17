# Paleo-detectors and neutrinos

*Is it possible find neutrinos from galactic core collapse SN with a paleo-detector?*

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb) [![DOI](https://zenodo.org/badge/142072044.svg)](https://zenodo.org/badge/latestdoi/142072044)  [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

Code for calculating track length spectra in paleo-detectors and exploring paleo-detector sensitivity to SN neutrinos. This code can be used to reproduce the results of [arXiv:1905.XXXX](http://arxiv.org/abs/1905.XXXX), "*Galactic Supernova Paleology*".

More information about paleo-detectors can also be found in [arXiv:1811.10549](http://arxiv.org/abs/1811.10549), [arXiv:1806.05991](http://arxiv.org/abs/1806.05991), and [arXiv:1811.06844](http://arxiv.org/abs/1811.06844).

**Authors:** Thomas D P Edwards, Bradley J Kavanagh, 

Please get in touch with any questions, comments or bug-reports.

### Overview: Paleopy

The core of the code runs on the [`paleopy`](https://github.com/tedwards2412/paleopy) package. This includes data for converting recoil energies to track lengths, along with tables of background distributions. This then allows you to calculate all the relevant track length distributions. The currently supported minerals are Nchwaningite, Sinjarite, Halite, Olivine, Gypsum and Phlogopite.

To run the notebooks it is required to install the paleopy package with:

    pip3 install git+https://github.com/tedwards2412/paleopy

Check out [`Notebooks/PlotSpectra.ipynb`](Notebooks/PlotSpectra.ipynb) for an illustration of how to use the code: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb)


### Notebooks

Details need to be added here


### Results

Tables of projected upper limits and discovery reach are output to [`ES/limits`](ES/limits).

### Requirements

The code in this repo should run with Python3. Standard library requirements are in [`requirements.txt`](requirements.txt). In addition, you will need [`swordfish`](https://github.com/cweniger/swordfish) for performing the statistical analysis and [`WIMpy`](https://github.com/bradkav/WIMpy_NREFT) for calculating the DM and neutrino spectra. 
