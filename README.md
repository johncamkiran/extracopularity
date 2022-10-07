-<img width="1080" alt="extracopularity-banner" src="https://user-images.githubusercontent.com/41493682/181186896-0c333843-a5bd-49c6-a6ff-51bd26bcb9ca.png">

# Extracopularity in 3D Particle Packings

Code for computing extracopularity coefficients in three-dimensional particle packings.

## Description

The [extracopularity coefficient](https://aip.scitation.org/doi/10.1063/5.0079985) is a local orientational order parameter for systems of interacting particles useful in the structural analysis of molecular dynamics simulations. In this repository, we provide implementations of our algorithm for computing extracopularity coefficients  in three-dimensional particle packings. A [MATLAB](https://www.mathworks.com/products/matlab.html) implementation is currently available, and a C++ implementation is under development. Our code accepts [LAMMPS](https://lammps.org/#gsc.tab=0) dump file inputs, which are easily visualized using [OVITO](https://www.ovito.org).

## Getting started

1. **Installation:** Place the function file *extracopularity.m* anywhere in the MATLAB [search path](https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html) or [present working directory](https://www.mathworks.com/help/matlab/ref/pwd.html).
2. **Execution:** Pass the full path of the target [LAMMPS dump file](https://docs.lammps.org/dump.html) as a string into the function. If the data is in another file format, simply [import](https://www.ovito.org/docs/current/usage/import.html#usage-import) it into OVITO and [export](https://www.ovito.org/manual/usage/export.html) it as a LAMMPS dump file before doing so.
3. **Visualization:** Import the output of the code into OVITO and use the [color coding modifier](https://www.ovito.org/docs/current/reference/pipelines/modifiers/color_coding.html) with the input property set to *c_extra*. Use the [slice modifier](https://www.ovito.org/docs/current/reference/pipelines/modifiers/slice.html#particles-modifiers-slice) to remove external surfaces if the simulation had been performed with periodic boundary conditions, as they are not recognized by current versions of the code.
4. **Postprocessing:** Use the [smooth trajectory modifier](https://www.ovito.org/docs/current/reference/pipelines/modifiers/smooth_trajectory.html#particles-modifiers-smooth-trajectory) in OVITO to reduce temporal noise. This must be placed lowest in the list of modifiers for it to produce the desired effect.

For further documentation, see the header comment of *extracopularity.m*

## Versions

The latest version performs "bond angle discretization" as per Ref. [2] if a "commonly encountered geometry" cannot be found as per Ref. [1]. A version that adheres strictly to Ref. [1] is preserved in the legacy branch. 

## Citation

If this software is used in the preparation of work that is to be published, please cite the following:

1. John Çamkıran, Fabian Parsch, and Glenn D. Hibbard , "A local orientational order parameter for systems of interacting particles", [J. Chem. Phys. 156, 091101 (2022)](https://doi.org/10.1063/5.0079985)

2. John Çamkıran, Fabian Parsch, and Glenn D. Hibbard , "On the topology of the space of coordination geometries", [arXiv:2207.12171 [math-ph] (2022)](https://arxiv.org/abs/2207.12171)

## Licenses

This software is distributed under the BSD-3-Clause License. For more details see the LICENSE file.
