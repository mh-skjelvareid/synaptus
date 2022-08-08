**Cite Synaptus** (DOI): [![DOI](https://joss.theoj.org/papers/10.21105/joss.04185/status.svg)](https://doi.org/10.21105/joss.04185)

**License:** [![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/mh-skjelvareid/synaptus/blob/master/LICENSE.md)


# Summary
Synaptus is a Matlab/Octave toolbox for synthetic aperture or array imaging. It was originally developed for ultrasonic imaging for non-destructive testing, but can be applied for similar imaging modes (e.g. ground penetrating radar). The toolbox focuses on algorithms implemented in the Fourier domain, and on imaging in multilayered structures (e.g. water, metal, rock).

The core functionality of the toolbox is to create focused images from raw (unfocused) pulse-echo data. Such data is produced by a transducer that transmits waves into a propagating medium, and records backscattered waves from within the medium. A backscattered "echo" is created when the waves interact with an object or layer with different physical properties than the propagating medium, e.g. a metal object in water. A measurement at a single point in space thus produces a 1-dimensional "depth profile". By moving the transducer laterally relative to the object under study, it is possible to create a 2- or 3-dimensional image of the object. A similar measurement can be performed using an array of multiple transducers. Due to the divergence of the transducer beams, the echoes from scattering objects are "smeared" laterally, making the images unfocused and hard to interpret. Examples of such images are given in the "Example raw and focused images" section below.

The algorithms in the Synaptus toolbox take such raw, unfocused images as input and processes them to create focused images. Conceptually, this is done by treating the raw images as measurements of a wave field, and using the wave equation to manipulate the wave field into a focused image. The details of this process is described in the PhD thesis included in the toolbox.


# Background
The toolbox was originally written by Martin H. Skjelvareid, as a collection of algorithms developed during his work as a PhD candidate. The PhD thesis is included in the "docs" folder, and can also be downloaded from https://hdl.handle.net/10037/4649 . As with so many things, "the devil is in the details" when it comes practical implementation of synthetic aperture algorithms. The toolbox is meant to help people who are new to the field and are looking to make implementations of published algorithms. The data sets included in the toolbox will hopefully be useful in the development of new algorithms for similar measurement geometries. The toolbox could also represent a collection of reference methods against which new algorithms are compared.

The name of the toolbox is an abbreviation of "Synthetic APerTure UltraSound".

# Example raw and focused images
The following images are taken from the thesis to illustrate some of the applications of the algorithms in the toolbox. The original figure numbering and captions have been included for context. See the `thesis.pdf` file under the 'docs' folder for further details.

## Bottom drilled holes in PMMA and aluminium blocks

<img src="graphics/AcrylicGlassAndAluminiumLayers_Setup.png" width="600">

<img src="graphics/AcrylicGlassAndAluminiumLayers_3DRender.png" width="550">

<img src="graphics/AcrylicGlassAndAluminiumLayers_RawData.png" width="550">

<img src="graphics/AcrylicGlassAndAluminiumLayers_Focused.png" width="650">


## Objects placed on cylindrical surface

<img src="graphics/ObjectsInPipe_Setup.png" width="550">

<img src="graphics/ObjectsInPipe_Images.png" width="600">


## Rusted pipe imaged from inside

<img src="graphics/RustedPipeImaging_Setup.png" width="600">

<img src="graphics/RustedPipeInterior.png" width="600">

<img src="graphics/RustedPipe_Slices_RawAndFocused.png" width="600">


# Algorithms
The main focus of the algorithms is on Fourier-domain synthetic aperture processing of ultrasound data. Fourier-domain processing is very common for synthetic aperture radar and sonar, but in the field of ultrasound, the time-domain "delay-and-sum" approach still dominates. With Fourier-domain processing, it is possible to extrapolate a sampled pulse-echo wavefield in both space and time ("wavefield migration"). One major advantage of this approach is that the wavefield is easily extrapolated between media with different wave velocities, enabling multi-layer imaging (very relevant for immersion ultrasound imaging).

The phase shift migration (PSM) algorithm works by migrating the recorded wavefield in small steps, and creating a focused image line/plane at each depth (using the "exploding reflector model"). From a processing standpoint, the method is not optimal, as the complete wavefield spectrum matrix has to be multiplied with a phase factor matrix at each step. However, since matrix multiplication is very fast in Matlab/Octave, the method is quite fast in practice.

The MULOK algorithm is a multi-layer version of the method known as Stolt migration, omega-k focusing, f-k focusing, etc. (several names exist for the same method). The wavefield extrapolation used in PSM is also used in MULOK, to extrapolate the wavefield between layers with different wave velocities. However, each layer is focused "all at once" by a resampling of the wavefield spectrum. This is generally a quite efficient method, but the practical performance very much depends on the efficiency of the interpolation method used (the linear interpolation mode of the interp1 function in Octave/Matlab is only moderately fast).

The CPSM algorithm is an adaptation of the PSM algorithm to a cylindrical imaging geometry (transducer pointing outward from cylindrical scanning surface). The solutions to the wave equation in a cylindrical geometry are Hankel functions. These can be used to extrapolate the wavefield, but are generally very time-consuming to compute. The algorithm also includes two alternative transfer functions which are approximate but much faster.



# Organization
The toolbox is organized into the following folders:
- 'core' contains the functions for synthetic aperture focusing. Each file represents a separate algorithm.
- 'datasets' contains datasets in .mat-format, used for test/demonstration of the algorithms
- 'test' contains test scripts for the algorithms
- 'misc' contains various function used to help in processing and plotting of results.
- 'learn' contains simplified versions of (some of) the algorithms in the toolbox, with additional plots of data at intermediate steps to help understanding.
- 'docs' contains relevant documentation (PhD thesis ++)
- 'experimental' contains "draft" code related to smaller concepts and ideas, including tilt compensation. The code is not fully polished/commented.

# Requirements
The toolbox requires a working base installation of Matlab or Octave, with some additional signal processing functions. Details for Matlab and Octave are given below.

## Matlab
Matlab is available for Linux, MacOS and Windows, and can be purchased from [Mathworks](https://mathworks.com/store/). Running Synaptus on Matlab requires:
- Base installation of Matlab (has been tested with version 2021a, installed on Windows)
- [Signal processing toolbox](https://se.mathworks.com/products/signal.html) (has been tested with version 8.6)

## Octave
GNU Octave is free and open source, and is available for both Linux, MacOS and Windows. See [GNU download and installation page](https://www.gnu.org/software/octave/download) for details.

Running Synaptus on Octave requires:
- Base installation of Octave (has been tested with version 5.2.0, installed on Ubuntu 20.04). On Debian / Ubuntu systems, Octave can be installed with `sudo apt-get install octave`.
- The [signal package](https://octave.sourceforge.io/signal/index.html) (has been tested with version 1.4.1). Octave installers for Windows bundle the signal package with the main application. On Debian/Ubuntu systems, the package can be installed using `sudo apt-get install octave-signal` Note that the signal package needs to be "loaded" inside Octave using the command `pkg load signal`. This can be automated by adding the command to a [Octave startup file](https://octave.org/doc/v6.4.0/Startup-Files.html#Startup-Files).


# Installation
Download the toolbox and add (at least) the "core" folder to the Matlab / Octave path. Run the scripts found under "test" to see example usage of the different algorithms. To run all the tests one after another, use the script "tests_runAll.m". Open and run the scripts under "learn" to see simplified versions of some of the algorithms, with plots.

# Documentation
The core algorithms are documented by function descriptions in the standard Matlab/Octave style. Use the `help` command to display documentation for a given function, e.g. `help psm`. The PhD thesis in the "docs" folder describes the theory behind the core algorithms.

The scripts in the "test" folder are meant to test the toolbox functionality, but also serve as an illustration of typical use of the toolbox on some example datasets.

The scripts in the "learn" folder contain simplified versions of the PSM and MULOK algorithms, with multiple plots showing the input data, intermediate steps, and the final focused image. The purpose of these scripts is to show the main steps of the algorithms, without the generalization/optimization "clutter" found in the core functions.

# Contributing
Contributions to the toolbox are most welcome; bug reports, suggestions for changes, datasets, new algorithms - anything you think is relevant. See [CONTRIBUTING](CONTRIBUTING.md) file for further details.

# License
[![License](https://img.shields.io/badge/License-GNU_GPLv3-orange.svg)](https://github.com/mh-skjelvareid/synaptus/blob/master/LICENSE.md)

# Citing Synaptus
If you use the Synaptus toolbox in your work, you should cite the following paper:

Skjelvareid, M. H., (2022). _Synaptus: A Matlab/Octave toolbox for synthetic aperture ultrasound imaging._ Journal of Open Source Software, 7(76), 4185, https://doi.org/10.21105/joss.04185

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04185/status.svg)](https://doi.org/10.21105/joss.04185)

# Acknowledgements
- The main elements of the toolbox were developed by M. H. Skjelvareid, working as an industrial PhD working at Breivoll Inspection Technologies (BIT) in Tromsø, Norway. The work was financed in equal parts by BIT and the Norwegian Research Council.
- M. H. Skjelvareid owes great thanks to Tomas Olofsson, who introduced him to phase shift migration and collaborated with him on multiple publications, and to Yngve Birkelund and Yngvar Larsen, who were his PhD advisors and co-authors.
