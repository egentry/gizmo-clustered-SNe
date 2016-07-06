#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * slightly (re-arranged, consolidated, and added compute_stellar_feedback and 
 * the gradients loop) by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void compute_grav_accelerations(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  if(ThisTask == 0)
    {
      printf("Start gravity force computation...\n");
      fflush(stdout);
    }


  gravity_tree();		/* computes gravity accel. */

  /* For the first timestep, we redo it to allow usage of 
   relative opening criterion for consistent accuracy */
  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();

  if(ThisTask == 0)
    {
      printf("gravity force computation done.\n");
      fflush(stdout);
    }
}



void compute_hydro_densities_and_forces(void)
{
  if(All.TotN_gas > 0)
    {
        if(ThisTask == 0)
        {
            printf("Start density & tree-update computation...\n"); fflush(stdout);
        }

        density();		/* computes density, and pressure */

#ifdef GRACKLE_FIX_TEMPERATURE
        if(RestartFlag == 0 && All.InitGasTemp > 0 && All.TimeStep == 0) fix_temperature();
#endif

#ifdef ADAPTIVE_GRAVSOFT_FORALL
        ags_density();
#endif
        force_update_hmax();	/* update kernel lengths in tree */
        /*! This function updates the hmax-values in tree nodes that hold SPH
         *  particles. These values are needed to find all neighbors in the
         *  hydro-force computation.  Since the Hsml-values are potentially changed
         *  in the SPH-denity computation, force_update_hmax() should be carried
         *  out before the hydrodynamical SPH forces are computed, i.e. after
         *  density().
         */
        
        if(ThisTask == 0)
        {
            printf("density & tree-update computation...\n"); fflush(stdout);
        }
        if(ThisTask == 0)
        {
            printf("Start gradient computation...\n"); fflush(stdout);
        }
        hydro_gradient_calc(); /* calculates the gradients of hydrodynamical quantities  */
        
        if(ThisTask == 0)
        {
            printf("gradient computation done.\n"); fflush(stdout);
        }
        
        if(ThisTask == 0)
        {
            printf("Start hydro-force computation...\n"); fflush(stdout);
        }
        
        hydro_force();		/* adds hydrodynamical accelerations and computes du/dt  */
        
        if(ThisTask == 0)
        {
            printf("hydro force computation done.\n"); fflush(stdout);
        }


    } else {
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        ags_density(); // if there are no gas particles but ags-all is active, still need to enter this loop //
        force_update_hmax();    /* update kernel lengths in tree */
#endif
    }
}


#ifdef GALSF
void compute_stellar_feedback(void)
{
    CPU_Step[CPU_MISC] += measure_time();

    /* first, check the mechanical sources of feedback */
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    mechanical_fb_calc(-1); /* compute weights for coupling */
    CPU_Step[CPU_SNIIHEATING] += measure_time();
#endif // (defined(FLAG_NOT_IN_PUBLIC_CODE)||defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE))
    
    /* alternatively use the 'turn off cooling' sub-grid feedback model */
    
    /* now do the local photo-ionization heating */
    
    /* finally (if we're not doing it in the star formation routine), do the local radiation pressure */
#if defined(GALSF_FB_RPWIND_FROMSTARS) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
    radiation_pressure_winds_consolidated();
    CPU_Step[CPU_LOCALWIND] += measure_time();
#endif
    
#ifdef GALSF_FB_LUPI
    //compute coupling for SNae and mass losses
    lupi_fb_calc(-1);
    //compute the actual SNa feedback
    if(All.FeedbackMode % 4 > 0) lupi_fb_calc(0);
    //compute IM mass star mass loss
    if(All.FeedbackMode > 4)     lupi_fb_calc(1);
#endif

#ifdef GENTRY_FB
    //compute the actual SNe feedback
    gentry_fb_calc();
    //compute IM mass star mass loss
    // if(All.FeedbackMode > 4)     lupi_fb_calc(1);
#endif
    CPU_Step[CPU_MISC] += measure_time();
}
#endif // GALSF //
