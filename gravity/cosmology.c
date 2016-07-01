#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"


/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO. The modifications
 * mostly center on added functionality for new modules, changed variables for 
 * cosmology, and consolidating the relevant functions into a single file.
 */


/* this function gets called regardless of the cosmology choices: 
    anything which modifies the growth history should live here. 
    This is the usual hubble function H0 * E(z); so for example 
    the proper time:
        dt = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
*/
double INLINE_FUNC hubble_function(double a)
{
    double hubble_a;
    
    hubble_a = All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a)
    + All.OmegaLambda;
    hubble_a = All.Hubble_H0_CodeUnits * sqrt(hubble_a);
    return (hubble_a);
}



/* Using Dark Energy instead of a Cosmological constant can be archived by
 * replacing Lambda by Lambda * a^(-3*(1+w)) in the Hubble function.
 * So easy to see that w = -1 gives back a standard Cosmological Constant !
 * Also w = -1/3 gives Lambda / a^2 which then cancel within the Hubble
 * function and is then equal to the dynamics of a universe with Lambda = 0 !
 *
 * For a time varying w once has to replace Lambda * a^(-3*(1+w)) by
 * Lambda * exp(Integral(a,1,3*(1+w)/a))
 *
 * One can now also read in colums for the change of the gravitational
 * constant and the correction by this to the hubble function.
 *
 * Additional once can read also an "external" hubble function from a
 * column of the dark energy file.
 *
 * Note that the first column is 'z+1' !
 *
 * Dark Energy does not alter the powerspectrum of initial conditions.
 * To get the same cluster for various values or functions of w, once
 * has do assign a new redshift to the initial cond. to match the
 * linear growth factors, so g(z=0)/g(z_ini) == g_w(z=0)/g_w(z_ini^new)
 * Also the initial velocities field has to be scaled by
 * (Hubble_w(z_ini^new)*Omega_w(z_ini^new)^0.6)/(Hubble(z_ini)*Omega(z_ini)^0.6)
 * where _w means the according functions including the terms for
 * Dark Energy.
 */


