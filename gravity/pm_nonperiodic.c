#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*! \file pm_nonperiodic.c
 *  \brief code for non-periodic FFT to compute long-range PM force
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * slightly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


#include "../allvars.h"
#include "../proto.h"

