# gizmo-clustered-SNe
3D simulations of clustered SNe with cooling, and with/without MHD
-------

Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   

Licensed under the GNU General Public License (v3)

-------

This work is an extension of my previous project, where I simulated some of the same star clusters, but with a 1D code that enforced spherical symmetry ([github.com/egentry/clustered_SNe](github.com/egentry/clustered_SNe)). For more details, you can check out the paper ([Gentry et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.465.2471G)). 

The general idea is to re-run some of those simulations for the same star clusters, but now in 3D, and maybe with magnetic fields.

-------

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
 - Any of the `runs/*/inputs/*test*` files. Those are qsub batch scripts.
 - Any of the `runs/*/inputs/*params.base` files. Files of those form are submitted as runtime parameter files for GIZMO. (In practice I don't actually submit `*params.base` -- I use them to build `*params.restart` files, adjusting things like `InitCondFile`, `TimeBegin`, `TimeMax`, `TimeOfFirstSnapshot` and `TimeBetSnapshot` for the appropriate values at a given part of the simulation.)
 - `scripts/generate_ICs.ipynb`. As the name suggests, it's the script that I use to generate the initial conditions for my simulations (uniform ISM, with/without magnetic fields, with varying box sizes and varying resolution).  At the bottom it also has a section "What to do next?" which will help you get going.

