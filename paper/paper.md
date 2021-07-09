---
title: 'Gala: A Python package for galactic dynamics'
tags:
  - ultrasound
  - synthetic aperture
  - SAFT
  - phased array
  - NDT
authors:
  - name: Martin H. Skjelvareid
    orcid: 0000-0003-2103-865X
    affiliation: 1
affiliations:
 - name: UiT - the Arctic University of Norway
   index: 1
date: 10 July 2021
bibliography: paper.bib

---

# Summary

Sonic imaging is performed by transmitting sound waves into an object and recording the waves scattered from within the object. In many applications, transmission and recording are performed with the same unit, a transducer. In some cases a single transducer is used, and by moving the transducer relative to the object of interest, a 2D or 3D map of reflections is created. Transducer arrays, which are basically several small transducers stacked along a line or on a grid, are also becoming increasingly common. By modulating the amplitude and time delay for transmission from each array element, arrays can create waves which are focused or shaped according to the application. The arrays also enable recording of scattered reflections at many different positions, without moving the array as a whole.

The recordings of scattered waves represent "raw data" for both single-transducer and array imaging setups. This data often suffers from poor resolution, making it hard (if not impossible) to interpret directly.  The data needs to be focused in order to create an image resembling the physical structure of an object. This process - creating focused images from raw pulse-echo ultrasound data - is the subject of the software presented here.

Sonic imaging is used in a number of different applications, such as geophysical exploration, medical imaging, and non-destructive testing (NDT) of industrial components. 



# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:

# Acknowledgements

(add acknowledgements)

# References