# gizmo-clustered-SNe
3D simulations of clustered SNe with cooling, and with/without MHD
-------

Current Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   
(Note: There have been many authors of GIZMO and GADGET before that. I'm only claiming to be the author of the modifications, setup and analysis that's specific to my project.)

Licensed under the GNU General Public License (v3)

-------

This work is an extension of my previous project, where I simulated some of the same star clusters, but with a 1D code that enforced spherical symmetry ([github.com/egentry/clustered_SNe](github.com/egentry/clustered_SNe)). For more details, you can check out the paper ([Gentry et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.465.2471G)). 

The general idea is to re-run some of those simulations for the same star clusters, but now in 3D, and maybe with magnetic fields.

-------
# Find your way here from the paper?
You might be interested in:
 - `scripts/publication_plots.ipynb`: shows how we created every plot in the paper. Includes png version of every plot. Also includes a bunch of extra plots / diagnostics related to the data being shown.
 - `scripts/plan_energy_content*.ipynb`: helps explain why we picked the times we did for our phase diagrams + snapshot density slices.  Also provides more information about the energy evolution of the simulations than what we could include in the paper.
 - `scripts/2d cooling curve.ipynb`: a 2d visualization of the cooling curve. Usually you just see the cooling curve at fixed density, but the density dependence does matter for some questions. This takes the same axes as our phase diagrams, and plots things like the specific cooling rate and the cooling time at each point in (density, temperature) space.  In theory, you can use take the mass-weighted phase diagram, multiply it by the 2d cooling curve, and get back the cooling-weighted phase diagram.


## Want data?
 - I can easily give you sqlite tables of reduced data for each simulation, containing things like mass, momentum, energy (kinetic, internal, magnetic) at each snapshot.
 - I can easily give you some versions of the raw 1D data.
 - The raw 3D data is _much_ more difficult. It is about 5 TB in size.  I might be able to run analysis scripts locally for you. Alternative we _might_ be able to work something out to transfer a subset of the data.  Email me if you want to talk about our options.

-------
# Want to actually run the code?

In general, the version of GIZMO that I'm using is only very slightly modified from the public version of GIZMO.  Most of the feedback injection that I'm doing actually happens outside of the C code---I stop the simulation at the time of each SN, create a new HDF5 snapshot file using python, then restart the C code.  This means the C code should be _fairly_ clean, although not perfectly so.  Get in touch with me over email if you have any questions. 

The next major change is that this version of GIZMO has Grackle cooling as an option. I got the initial implementation from Alessandro Lupi, but I've since made it more standardized (rather than his custom version of Grackle) and I've upgraded it to v3 of the Grackle API
. So in order to use Grackle, you should install it as a stand-alone library, then link it as normal.

The final change to the c code is that I've disabled the creation of restartfiles, because they take up a ton of memory and tended to break my code more often than they were useful. If you'd like to revert the c code to create restartfiles, then you should remove the `return;` line in the `restart` function in the `restart.c` file.)

For more details about the project see:
 - Phil Hopkin's [GIZMO page](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html)
 - the README for the public version of GIZMO (saved as `README-Hopkins.md` in this repo)
 - the [Grackle documentation](grackle.readthedocs.io)

I've also added a bunch of python scripts/jupyter notebooks to help generate initial conditions, add SN feedback particles, and visualize the results.  If you're curious, you can check them out in the `scripts` directory.

If you're wondering where to start, here are a few key files to look at:
 - `scripts/generate_ICs.ipynb`. As the name suggests, it's the script that I use to generate the initial conditions for my simulations (uniform ISM, with/without magnetic fields, with varying box sizes and varying resolution).  At the bottom it also has a section "What to do next?" which will help you get going.
 - Any of the `runs/*/inputs/*_loop-*` files. Those are qsub batch scripts.
 - Any of the `runs/*/inputs/*params.base` files. Files of this form are submitted as input runtime parameter files for GIZMO. (In practice I don't actually submit `*params.base` -- I use them to build `*params.restart` files, adjusting things like `InitCondFile`, `TimeBegin`, `TimeMax`, `TimeOfFirstSnapshot` and `TimeBetSnapshot` for the appropriate values at a given part of the simulation.)

# Requirements
## For the core C code
More explanation [here](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html#tutorial-requirements)
 - MPI-enabled C compiler
 - GSL
 - HDF5
 - FFTW 2.x (**not** 3.x); MPI-capable
 - Grackle (version 3)
 
 ## For the python code
 This will only be for the notebooks to generate initial conditions and for the scripts to add SNe. This will not include all the required packages for the various analysis notebooks.
  - Python 3
  - numpy
  - astropy
  - h5py
  - scipy
