#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"

#include "./cooling.h"


/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM, 
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including 
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these 
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel (volker.springel@h-its.org). The code has been modified heavily by 
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; everything except the original metal-free free-free and 
 *   photo-ionization heating physics has been added (or re-written), and the iteration routine to converge to 
 *   temperatures has been significantly modified.
 */


#ifdef COOLING


static double XH = HYDROGEN_MASSFRAC;	/* hydrogen abundance by mass */

#define eV_to_K   11606.0
#define eV_to_erg 1.60184e-12


/* this is just a simple loop if all we're doing is cooling (no star formation) */
void cooling_only(void)
{
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            do_the_cooling_for_particle(i);
        } // if(P[i].Type == 0 && P[i].Mass > 0)
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
} // void cooling_only(void)





/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle(int i)
{
    double unew;
    double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
    double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
    if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
        
        double ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */
        double uold = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);
        
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        double grav_acc; int k;
        for(k = 0; k < 3; k++)
        {
            grav_acc = All.cf_a2inv * P[i].GravAccel[k];
            SphP[i].DtInternalEnergy -= SphP[i].GravWorkTerm[k] * All.cf_atime * grav_acc;
        }
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
        double du = SphP[i].DtInternalEnergy * dtime;
        if(du < -0.5*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy = -0.5*SphP[i].InternalEnergy / dtime;}
        if(du >  50.*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy =  50.*SphP[i].InternalEnergy / dtime;}
        /* and convert to cgs before use in the cooling sub-routine */
        SphP[i].DtInternalEnergy *= All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH);
        
        
        /* Call the actual COOLING subroutine! */
        unew = DoCooling(uold, SphP[i].Density * All.cf_a3inv, dtime, &ne, i);
        
        
        
        
        

        
        
        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if FLAG_NOT_IN_PUBLIC_CODE
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated. 
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        SphP[i].InternalEnergy = unew;
        SphP[i].Ne = ne;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
        SphP[i].DtInternalEnergy = 0;
        
        
        
    } // closes if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}






/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling(double u_old, double rho, double dt, double *ne_guess, int target)
{
  double u, du;

#ifdef GRACKLE
#if !defined(FLAG_NOT_IN_PUBLIC_CODE) && !defined(GRACKLE_FULLYIMPLICIT)
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below
     to split the step */
    du = dt * SphP[target].DtInternalEnergy / (All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH));
    u_old += 0.5*du;
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
    /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
     term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
     the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
     to which we then add the heating term from hydro and compute the solution over a full timestep */
    double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(r-1)/log(r);} r=DMAX(0,DMIN(r,1));
    du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
#if (GRACKLE_CHEMISTRY>=2)
    SphP[target].Gamma = CallGrackle(u, rho, 0, ne_guess, target, 4);
#endif
    return DMAX(u,All.MinEgySpec);
#endif
    
    
}



/* returns cooling time. 
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess, int target)
{
    double u;
    double LambdaNet, coolingtime;

#if defined(GRACKLE) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
    coolingtime = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(coolingtime >= 0) coolingtime = 0.0;
    coolingtime *= All.HubbleParam / All.UnitTime_in_s;
    return coolingtime;
#endif
    
}




/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess, int target)
{
    
}

/*  this function computes the self-consistent temperature and electron fraction */ 
double ThermalProperties(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer, double *mu_pointer, int target)
{
 return 0;
}




void InitCoolMemory(void)
{

}



void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19 
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
}





/* table input (from file TREECOOL) for ionizing parameters */
/* NOTE: we've switched to using the updated TREECOOL from CAFG, june11 version */


void ReadIonizeParams(char *fname)
{
}


void IonizeParams(void)
{
  IonizeParamsTable();

  /*
     IonizeParamsFunction();
   */
}



void IonizeParamsTable(void)
{
  return;
}


void SetZeroIonization(void)
{
}


void IonizeParamsFunction(void)
{

}



void InitCool(void)
{
    if(ThisTask == 0)
        printf("Initializing cooling ...\n");
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();
    
#ifdef GRACKLE
    InitGrackle();
#endif
    
    InitCoolMemory();
    MakeCoolingTable();
    ReadIonizeParams("TREECOOL");
    IonizeParams();
}












#endif
