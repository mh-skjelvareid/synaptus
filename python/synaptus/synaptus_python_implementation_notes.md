# Synaptus-python implementation notes

## Frequency sign conventions
The original version of synaptus was based on using negative omega values. This now
seems like an unnecessary complication, as wavefield extrapolation can be done using
positive omega as well. Conventions differ, but it seems that algorithms related to
synthetic aperture imaging and ultrasound usually use positive omega _and_ the positive
root when computing k_z, and then apply a negative sign in the final phase shift.

Note that approach is already used in the array_psm algorithm.

## Inheritance structure
Some attributes and methods are the same for PSM, MULOK, and CPSM. What is the best way
to avoid code duplication? Using a base class (abstract or regular) is an obvious
solution - but what should it include?

Common attributes:
- raw_data
- ndim, nt, nx, ny
- fs, f_low, f_high
- t_delay
- x_step, y_step (only for PSM and MULOK)
- sound_velocities
- layer_thicknesses
- nfft_t 
- nfft_x, nfft_y (only for PSM and MULOK)
- omega_vec 
- kx_vec, ky_vec (only for PSM and MULOK)
- OMEGA, KX, KZ

The cylindrical coordinate system of CPSM is a bit of a headache in that the _names_ of
all the terms used are different (x,y,z vs r, phi, z, and Fourier-domain counterparts).
The role of the z axis is also different.  

The CPSM algorithms is kind of av niche application - focusing on PSM and MULOK makes
sense. However, the elements that _are_ the same for PSM, MULOK and CPSM can be made
independent functions, placed outside the class definitions (in utils?). 

### Independent functions
- nextpow2
- get_time_vec (incl. t_delay)
- get_omega_vec
- get_passband_indices
- get_grids


### Base class: MultilayerCartesianPulseEchoDataset
Class containing a dataset with acquisition parameters, time-space coordinate axes, a
Fourier-domain wavefield dataset, and corresponding omega-k coordinate axes. 

Methods:
- get_omega_k_grids
- time-shift wavefield
- extrapolate_wavefield()

### PSM
Inherits from MultilayerCartesianPulseEchoDataset

Methods:
- get_phase_shift_tensor
- phase_shift_migrate

### MULOK
Inherits from MultilayerCartesianPulseEchoDataset

Methods:
- interpolate_wavefield



