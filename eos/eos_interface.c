#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gsl/gsl_const_cgsm.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "eos.h"
#include "../allvars.h"


#define BITMASK_SET_FLAG(BITMASK,FLAG)      (BITMASK) |= (FLAG)
#define BITMASK_SET_ALL_FLAGS(BITMASK)      (BITMASK = ~(0))
#define BITMASK_UNSET_FLAG(BITMASK,FLAG)    (BITMASK) &= ~(FLAG)
#define BITMASK_UNSET_ALL_FLAGS(BITMASK)    (BITMASK) = 0
#define BITMASK_CHECK_FLAG(BITMASK,FLAG)    (((BITMASK) & (FLAG)) == (FLAG))

#define EOS_ERR_VALID         0
#define EOS_ERR_RHO_LT_RHOMIN 1
#define EOS_ERR_RHO_GT_RHOMAX 2
#define EOS_ERR_EPS_LT_EPSMIN 4
#define EOS_ERR_EPS_GT_EPSMAX 8
#define EOS_ERR_COMPOSITION   16

static int eos_input_to_cgs(struct eos_input * vars);
static int eos_output_from_cgs(struct eos_output * vars);

static int eos_validate(struct eos_input const * vars, struct eos_input * vars_adj, int * bitmask);
static int eos_compute_from_valid(struct eos_input const * in, struct eos_output * out);

#ifdef EOS_TABULATED
int eos_init(char const * eos_table_fname)
{
  if(access(eos_table_fname, R_OK) != 0)
  {
    fprintf(stderr, "Could not read \"%s\"\n", eos_table_fname);
    return 1;
  }
  return 0;
}

int eos_cleanup()
{
  return 0;
}
#endif

int eos_compute(struct eos_input const * in_, struct eos_output * out_)
{
  struct eos_input in, in_adj;
  memcpy(&in, in_, sizeof(in));
#ifdef EOS_USES_CGS
  eos_input_to_cgs(&in);
#endif

  int bitmask = 0;
  int ierr = eos_validate(&in, &in_adj, &bitmask);
  assert(!ierr);
  if(bitmask != EOS_ERR_VALID)
  {
    fprintf(stderr, "EOS ERROR:");
    if(BITMASK_CHECK_FLAG(bitmask, EOS_ERR_COMPOSITION))
      fprintf(stderr, "/invalid composition");
    if(BITMASK_CHECK_FLAG(bitmask, EOS_ERR_RHO_LT_RHOMIN))
      fprintf(stderr, "/density too low");
    if(BITMASK_CHECK_FLAG(bitmask, EOS_ERR_RHO_GT_RHOMAX))
      fprintf(stderr, "/density too large");
    if(BITMASK_CHECK_FLAG(bitmask, EOS_ERR_EPS_LT_EPSMIN))
      fprintf(stderr, "/temperature too low");
    if(BITMASK_CHECK_FLAG(bitmask, EOS_ERR_EPS_GT_EPSMAX))
      fprintf(stderr, "/temperature too high");
    fprintf(stderr, "\n");

#ifdef EOS_USES_CGS
    char const * unit_dens = "g/cm^3";
    char const * unit_ene  = "erg/g";
#else
    char const * unit_dens = "";
    char const * unit_ene  = "";
#endif
     
    fprintf(stderr, "  rho  = %.19e %s\n", in.rho, unit_dens);
    fprintf(stderr, "  eps  = %.19e %s\n", in.eps, unit_ene);
#ifdef EOS_CARRIES_YE
    fprintf(stderr, "  Ye   = %.19e\n", in.Ye);
#endif
#ifdef EOS_CARRIES_ABAR
    fprintf(stderr, "  Abar = %.19e\n", in.Abar);
#endif
    fprintf(stderr, "Using 0th order extrapolation\n");
    memcpy(&in, &in_adj, sizeof(in));
  }
  struct eos_output out;
  ierr = eos_compute_from_valid(&in, &out);
  assert(!ierr);
#ifdef EOS_USES_CGS
  ierr = eos_output_from_cgs(&out);
  assert(!ierr);
#endif
  memcpy(out_, &out, sizeof(out));

  return 0;
}

static int eos_input_to_cgs(struct eos_input * vars)
{
  vars->rho *= All.UnitDensity_in_cgs;
  vars->eps *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;
  return 0;
}

static int eos_output_from_cgs(struct eos_output * vars)
{
  vars->press  /= All.UnitPressure_in_cgs;
  vars->csound /= All.UnitVelocity_in_cm_per_s;
  return 0;
}

static int eos_validate(struct eos_input const * vars, struct eos_input * vars_adj, int * bitmask)
{
  *bitmask = EOS_ERR_VALID;
  memcpy(vars_adj, vars, sizeof(*vars));


  return 0;
}

static int eos_compute_from_valid(struct eos_input const * in, struct eos_output * out)
{

  return 0;
}
