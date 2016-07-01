#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"

#ifdef GRACKLE
#include <grackle.h>
#define ENDRUNVAL 91234

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return pressure
//     4 == calculate and return gamma (only valid when GRACKLE_CHEMISTRY>0)
//
double CallGrackle(double u_old, double rho, double dt, double *ne_guess, int target, int mode)
{
    gr_float returnval = 0.0;
    gr_float edot = 0.0;
#if defined(GRACKLE_FULLYIMPLICIT) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
    edot = SphP[target].DtInternalEnergy / (All.HubbleParam * All.UnitEnergy_in_cgs /
		(All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/HYDROGEN_MASSFRAC)) * rho;
#endif
    
    // Set grid dimension and size.
    // grid_start and grid_end are used to ignore ghost zones.
    int field_size = 1;
    int grid_rank = 3;
    int grid_dimension[3], grid_start[3], grid_end[3];
    int i;
    for (i = 0;i < 3;i++) {
        grid_dimension[i] = 1; // the active dimension not including ghost zones.
        grid_start[i]     = 0;
        grid_end[i]       = 0;
    }
    grid_dimension[0] = field_size;
    grid_end[0]       = field_size - 1;
    
    gr_float density, metal_density, energy, velx, vely, velz;
    gr_float cooling_time, temperature, pressure, gamma;
    velx          = SphP[target].VelPred[0];
    vely          = SphP[target].VelPred[1];
    velz          = SphP[target].VelPred[2];
    density       = rho;
    energy        = u_old;
#ifdef GRACKLE_OPTS
    metal_density = density * P[target].Metallicity[0];
#else
    metal_density = density * 0.02;
#endif
    gamma         = GAMMA;
    
    gr_float ne_density;
    gr_float HI_density, HII_density, HM_density;
    gr_float HeI_density, HeII_density, HeIII_density;
    gr_float H2I_density, H2II_density;
    gr_float DI_density, DII_density, HDI_density;
    gr_float tiny = 1.0e-20;
    
    ne_density    = density * tiny;
    
    HI_density    = density * tiny;
    HII_density   = density * tiny;
    HM_density    = density * tiny;
    
    HeI_density   = density * tiny;
    HeII_density  = density * tiny;
    HeIII_density = density * tiny;

    H2I_density   = density * tiny;
    H2II_density  = density * tiny;
    DI_density    = density * tiny;
    DII_density   = density * tiny;
    HDI_density   = density * tiny;
    
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
    // Atomic
    ne_density    = density * *ne_guess;
    
    HI_density    = density * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
    HII_density   = density * SphP[target].grHII;
    HM_density    = density * SphP[target].grHM;
    
    HeI_density   = density * SphP[target].grHeI;
    HeII_density  = density * SphP[target].grHeII;
    HeIII_density = density * SphP[target].grHeIII;
#endif
    
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    H2I_density  = density * SphP[target].grH2I;
    H2II_density = density * SphP[target].grH2II;
    gamma = SphP[target].Gamma;
#endif
    
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    DI_density   = density * SphP[target].grDI;
    DII_density  = density * SphP[target].grDII;
    HDI_density  = density * SphP[target].grHDI;
#endif
    
    switch(mode) {
        case 0:  //solve chemistry & update values
            if(solve_chemistry(&All.GrackleUnits,
                               All.cf_atime, dt,
                               grid_rank, grid_dimension,
                               grid_start, grid_end,
                               &density, &energy,
                               &velx, &vely, &velz,
                               &HI_density, &HII_density, &HM_density,
                               &HeI_density, &HeII_density, &HeIII_density,
                               &H2I_density, &H2II_density,
                               &DI_density, &DII_density, &HDI_density,
                               &ne_density, &metal_density,edot) == 0) {
                fprintf(stderr, "Error in solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
            // Assign variables back
            *ne_guess            = ne_density    / density;
            
            SphP[target].grHI    = HI_density    / density;
            SphP[target].grHII   = HII_density   / density;
            SphP[target].grHM    = HM_density    / density;
            
            SphP[target].grHeI   = HeI_density   / density;
            SphP[target].grHeII  = HeII_density  / density;
            SphP[target].grHeIII = HeIII_density / density;
#endif
            
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grH2I   = H2I_density   / density;
            SphP[target].grH2II  = H2II_density  / density;
#endif
            
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = DI_density    / density;
            SphP[target].grDII   = DII_density   / density;
            SphP[target].grHDI   = HDI_density   / density;
#endif
            returnval = energy;
            break;
            
        case 1:  //cooling time
            if(calculate_cooling_time(&All.GrackleUnits, All.cf_atime,
                                      grid_rank, grid_dimension,
                                      grid_start, grid_end,
                                      &density, &energy,
                                      &velx, &vely, &velz,
                                      &HI_density, &HII_density, &HM_density,
                                      &HeI_density, &HeII_density, &HeIII_density,
                                      &H2I_density, &H2II_density,
                                      &DI_density, &DII_density, &HDI_density,
                                      &ne_density, &metal_density, edot,
                                      &cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            if(calculate_temperature(&All.GrackleUnits, All.cf_atime,
                                     grid_rank, grid_dimension,
                                     grid_start, grid_end,
                                     &density, &energy,
                                     &HI_density, &HII_density, &HM_density,
                                     &HeI_density, &HeII_density, &HeIII_density,
                                     &H2I_density, &H2II_density,
                                     &DI_density, &DII_density, &HDI_density,
                                     &ne_density, &metal_density,
                                     &temperature) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure
            if(calculate_pressure(&All.GrackleUnits, All.cf_atime,
                                  grid_rank, grid_dimension,
                                  grid_start, grid_end,
                                  &density, &energy,
                                  &HI_density, &HII_density, &HM_density,
                                  &HeI_density, &HeII_density, &HeIII_density,
                                  &H2I_density, &H2II_density,
                                  &DI_density, &DII_density, &HDI_density,
                                  &ne_density, &metal_density,
                                  &pressure) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
        case 4:  //calculate gamma
            if(calculate_gamma(&All.GrackleUnits, All.cf_atime,
                               grid_rank, grid_dimension,
                               grid_start, grid_end,
                               &density, &energy,
                               &HI_density, &HII_density, &HM_density,
                               &HeI_density, &HeII_density, &HeIII_density,
                               &H2I_density, &H2II_density,
                               &DI_density, &DII_density, &HDI_density,
                               &ne_density, &metal_density,
                               &gamma) == 0) {
                fprintf(stderr, "Error in calculate_gamma.\n");
                endrun(ENDRUNVAL);
            }
            returnval = gamma;
            break;
    } //end switch
    
    return returnval;
}




//Initialize Grackle
void InitGrackle(void)
{
    grackle_verbose = 0;
    // Enable output
    if(ThisTask == 0) grackle_verbose = 1;
    
    // First, set up the units system.
    // These are conversions from code units to cgs.
    // The following parameter must be kept equal to 0 ALWAYS!!!
    All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm / All.HubbleParam;
    All.GrackleUnits.time_units           = All.UnitTime_in_s / All.HubbleParam;
    All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
    All.GrackleUnits.a_units              = 1.0; // units for the expansion factor
    
    // Second, create a chemistry object for parameters and rate data.
    if (set_default_chemistry_parameters() == 0) {
        fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
        exit(ENDRUNVAL);
    }
    // Third, set parameter values for chemistry & cooling
    
    /* optional flags:: */
    
    // Flag to control which three-body H2 formation rate is used.
    //    0: Abel, Bryan & Norman (2002),
    //    1: Palla, Salpeter & Stahler (1983),
    //    2: Cohen & Westberg (1983),
    //    3: Flower & Harris (2007),
    //    4: Glover (2008).
    //    These are discussed in Turk et. al. (2011). Default: 0.
    grackle_data.three_body_rate        = 0;
    
#ifdef GRACKLE_OPTS
    grackle_data.metal_cooling          = All.MetalCooling;    // metal cooling on
#else
    grackle_data.metal_cooling          = 0;                   // metal cooling on
#endif
    grackle_data.h2_on_dust             = 0;                   // dust cooling/chemistry off
    grackle_data.photoelectric_heating            = 0;
    grackle_data.photoelectric_heating_rate       = 8.5e-26;
    
    // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    grackle_data.cmb_temperature_floor  = 1;
    // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
#ifdef GRACKLE_OPTS
    grackle_data.UVbackground           = All.UVBackgroundOn;   // UV background on
#else
    grackle_data.UVbackground           = 0;                  // UV background on
#endif
    // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    grackle_data.Compton_xray_heating   = 1;
    
    
    // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    grackle_data.cie_cooling                      = 0;
    // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    grackle_data.h2_optical_depth_approximation   = 0;
    
    // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
    //    in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    grackle_data.LWbackground_intensity           = 0;
    // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
    //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    grackle_data.LWbackground_sawtooth_suppression = 0;
    
    
    /* fixed flags:: */
    
    // Flag to activate the grackle machinery:
    grackle_data.use_grackle            = 1;                   // grackle on (duh)
    // Path to the data file containing the metal cooling and UV background tables:
    grackle_data.grackle_data_file      = All.GrackleDataFile; // data file
    // Flag to include radiative cooling and actually update the thermal energy during the
    // chemistry solver. If off, the chemistry species will still be updated. The most
    // common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    grackle_data.with_radiative_cooling = 1;                   // cooling on
    // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    grackle_data.Gamma                  = GAMMA;              // our eos set in Config.sh
    // Flag to control which primordial chemistry network is used (set by Config file)
#ifndef GRACKLE_CHEMISTRY
    grackle_data.primordial_chemistry = 0;                     // fully tabulated cooling
#else
    grackle_data.primordial_chemistry = GRACKLE_CHEMISTRY;
#endif
    
    // Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    double a_value = 1.0;
    if(All.ComovingIntegrationOn) a_value = All.TimeBegin;
    
    // Finally, initialize the chemistry object.
    if (initialize_chemistry_data(&All.GrackleUnits, a_value) == 0) {
        fprintf(stderr, "Error in initialize_chemistry_data.\n");
        exit(ENDRUNVAL);
    }
    
    if(ThisTask == 0)
        printf("GRACKLE INITIALIZED\n");
}

#ifdef GRACKLE_FIX_TEMPERATURE
void fix_temperature()
{
        double temp_units=1./(BOLTZMANN/PROTONMASS/(GAMMA-1)/pow(All.UnitVelocity_in_cm_per_s,2.0));
        double temp;
        double mu;
        int i;
        for(i=0;i<NumPart;i++)
        {
          if(P[i].Type!=0) continue;
          temp=All.InitGasTemp;
          mu=0.59;
          SphP[i].InternalEnergy=temp/temp_units/mu;
          do{
            SphP[i].InternalEnergy=All.InitGasTemp/temp*SphP[i].InternalEnergy;

            temp=CallGrackle(SphP[i].InternalEnergy,SphP[i].Density,0,&(SphP[i].Ne),i,2);
            mu=temp/SphP[i].InternalEnergy/temp_units;
          }while(fabs(temp-All.InitGasTemp)>1.0e-2);
          SphP[i].InternalEnergy=temp/temp_units/mu;
          SphP[i].InternalEnergyPred=temp/temp_units/mu;
	  SphP[i].Pressure = get_pressure(i);

        }
}
#endif

#if (GRACKLE_CHEMISTRY >= 1)
void init_species()
{
  //if(!All.CoolingOn) return;
  FILE* fspecies=fopen("PrimChemData.dat","r");
  double species[12];
  char line[200];
  char* spname;
  int i,j;
  if(fspecies!=NULL)
  {
    if(ThisTask==0) printf("Species abundancies from data file\n");
    for(i=0;i<12;i++)
    {
        spname=fgets(line,80,fspecies);
        spname=strtok(line,"=");
        spname=strtok(NULL,"");
        species[i]=atof(spname);
    }
    fclose(fspecies);
  }
  else
  {
        if(ThisTask==0) printf("No PrimChemData.dat file found! Species will be initialized in \
neutral gas conditions\n");
        fflush(stdout);
    for(i=0;i<12;i++)
        species[i]=1.0e-20;
    species[0]=HYDROGEN_MASSFRAC;//grackle_data.HydrogenFractionByMass;
    species[2]=(1-HYDROGEN_MASSFRAC);//grackle_data.HydrogenFractionByMass);
    //species[9]=grackle_data.DeuteriumToHydrogenRatio*species[0];
  }

  if(ThisTask==0)
  {
        printf("-------------------------\n");
        printf("- HI  %5.3e / HII  %5.3e / Ne    %5.3e -\n",species[0],species[1],species[5]);
        printf("- HeI %5.3e / HeII %5.3e / HeIII %5.3e -\n",species[2],species[3],species[4]);
        printf("- HM  %5.3e / H2I  %5.3e / H2II  %5.3e -\n",species[6],species[7],species[8]);
        printf("- DI  %5.3e / DII  %5.3e / HDI   %5.3e -\n",species[9],species[10],species[11]);
  }
  fflush(stdout);

  for(i=0;i<NumPart;i++)
  {
       if(P[i].Type==0)
       {
             SphP[i].Ne      =DMAX(species[5] ,1.0e-20);
             SphP[i].grHI   =DMAX(species[0] ,1.0e-20);
             SphP[i].grHII  =DMAX(species[1] ,1.0e-20);
             SphP[i].grHeI  =DMAX(species[2] ,1.0e-20);
             SphP[i].grHeII =DMAX(species[3] ,1.0e-20);
             SphP[i].grHeIII=DMAX(species[4] ,1.0e-20);
             SphP[i].grHM   =DMAX(species[6] ,1.0e-20);
#if (GRACKLE_CHEMISTRY>=2)
             SphP[i].grH2I  =DMAX(species[7] ,1.0e-20);
             SphP[i].grH2II =DMAX(species[8] ,1.0e-20);
   	     SphP[i].Gamma=GAMMA;
#endif
#if (GRACKLE_CHEMISTRY>=3)
             SphP[i].grDI   =DMAX(species[9] ,1.0e-20);
             SphP[i].grDII  =DMAX(species[10],1.0e-20);
             SphP[i].grHDI  =DMAX(species[11],1.0e-20);
#endif
       }
  }
}
#endif


#endif  //GRACKLE
