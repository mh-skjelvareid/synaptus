## Synaptus
Synaptus is a Python and Matlab/Octave toolbox for synthetic aperture and array imaging.
It was originally developed for ultrasonic imaging for non-destructive testing, but can
be applied for similar imaging modes (e.g. ground penetrating radar). The toolbox
focuses on algorithms implemented in the Fourier domain, and on imaging in multilayered
structures (e.g. water, metal, rock).

The core functionality of the toolbox is to create focused images from raw (unfocused)
pulse-echo data. Such data is produced by a transducer that transmits waves into a
propagating medium, and records backscattered waves from within the medium. A
backscattered "echo" is created when the waves interact with an object or layer with
different physical properties than the propagating medium, e.g. a metal object in water.
A measurement at a single point in space thus produces a 1-dimensional "depth profile".
By moving the transducer laterally relative to the object under study, it is possible to
create a 2- or 3-dimensional image of the object. A similar measurement can be performed
using an array of multiple transducers. Due to the divergence of the transducer beams,
the echoes from scattering objects are "smeared" laterally, making the images unfocused
and hard to interpret. Examples of such images are given in the "Example raw and focused
images" section below.

The algorithms in the Synaptus toolbox take such raw, unfocused images as input and
processes them to create focused images. Conceptually, this is done by treating the raw
images as measurements of a wave field, and using the wave equation to manipulate the
wave field into a focused image. The details of this process is described in the PhD
thesis included in the toolbox.
