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

Beyond that, the only non-standard part of the C code is the inclusion of Grackle (v2) to calculate metallicity-dependent cooling.  It its current form (circa 6 Nov 2017) it assumes a custom version of Grackle (although I'm not actually using any of the customizations), which has been included within this repo.  Maybe in the future I'll tear out those customizations, and upgrade to the Grackle v3 interface.  But, if you're reading this, then I guess I haven't gotten around to it yet.

For more details about the project see:
 - Phil Hopkin's [GIZMO page](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html)
 - the README for the public version of GIZMO (saved as `README-Hopkins.md` in this repo)
 - the [Grackle documentation](grackle.readthedocs.io)

I've also added a bunch of python scripts/jupyter notebooks to help generate initial conditions, add SN feedback particles, and visualize the results.  If you're curious, you can check them out in the `scripts` directory.
