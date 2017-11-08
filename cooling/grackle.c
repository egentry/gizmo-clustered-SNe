#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <wordexp.h>
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

//     gr_float edot = 0.0;
// #if defined(GRACKLE_FULLYIMPLICIT) && !defined(FLAG_NOT_IN_PUBLIC_CODE)
//     edot = SphP[target].DtInternalEnergy / (All.HubbleParam * All.UnitEnergy_in_cgs /
// 		(All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/HYDROGEN_MASSFRAC)) * rho;
// #endif

    gr_float cooling_time, temperature, pressure, gamma;

    All.GrackleUnits.a_value = All.cf_atime;

    grackle_field_data my_fields;
    
    // Set grid dimension and size.
    // grid_start and grid_end are used to ignore ghost zones.
    const int field_size = 1;
    my_fields.grid_rank = 3;
    my_fields.grid_dimension = malloc(field_size * sizeof(int));
    my_fields.grid_start     = malloc(field_size * sizeof(int));
    my_fields.grid_end       = malloc(field_size * sizeof(int));

    int i;
    for (i = 0;i < 3;i++) {
      my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
      my_fields.grid_start[i] = 0;
      my_fields.grid_end[i] = 0;
    }

    my_fields.grid_dimension[0] = field_size;
    my_fields.grid_end[0] = field_size - 1;

    my_fields.density         = malloc(field_size * sizeof(gr_float));
    my_fields.internal_energy = malloc(field_size * sizeof(gr_float));
    my_fields.x_velocity      = malloc(field_size * sizeof(gr_float));
    my_fields.y_velocity      = malloc(field_size * sizeof(gr_float));
    my_fields.z_velocity      = malloc(field_size * sizeof(gr_float));

    my_fields.x_velocity[0]      = SphP[target].VelPred[0];
    my_fields.y_velocity[0]      = SphP[target].VelPred[1];
    my_fields.z_velocity[0]      = SphP[target].VelPred[2];
    my_fields.density[0]         = rho;
    my_fields.internal_energy[0] = u_old;

    // for metal_cooling = 1
    my_fields.metal_density    = malloc(field_size * sizeof(gr_float));
#ifdef GRACKLE_OPTS
    my_fields.metal_density[0] = rho * P[target].Metallicity[0];
#else
    my_fields.metal_density[0] = rho * 0.02;
#endif
    gamma         = GAMMA;

    // for primordial_chemistry >= 1
    my_fields.HI_density      = malloc(field_size * sizeof(gr_float));
    my_fields.HII_density     = malloc(field_size * sizeof(gr_float));
    my_fields.HeI_density     = malloc(field_size * sizeof(gr_float));
    my_fields.HeII_density    = malloc(field_size * sizeof(gr_float));
    my_fields.HeIII_density   = malloc(field_size * sizeof(gr_float));
    my_fields.e_density       = malloc(field_size * sizeof(gr_float));
    // for primordial_chemistry >= 2
    my_fields.HM_density      = malloc(field_size * sizeof(gr_float));
    my_fields.H2I_density     = malloc(field_size * sizeof(gr_float));
    my_fields.H2II_density    = malloc(field_size * sizeof(gr_float));
    // for primordial_chemistry >= 3
    my_fields.DI_density      = malloc(field_size * sizeof(gr_float));
    my_fields.DII_density     = malloc(field_size * sizeof(gr_float));
    my_fields.HDI_density     = malloc(field_size * sizeof(gr_float));
    
    const double tiny_number = 1.0e-20;

    // bunch of default values, to be overwritten as needed
    my_fields.HI_density[0]    = grackle_data->HydrogenFractionByMass * my_fields.density[0];
    my_fields.HII_density[0]   = tiny_number  * my_fields.density[0];
    my_fields.HM_density[0]    = tiny_number  * my_fields.density[0];
    my_fields.HeI_density[0]   = (1.0 - grackle_data->HydrogenFractionByMass) * my_fields.density[0];
    my_fields.HeII_density[0]  = tiny_number  * my_fields.density[0];
    my_fields.HeIII_density[0] = tiny_number  * my_fields.density[0];
    my_fields.H2I_density[0]   = tiny_number  * my_fields.density[0];
    my_fields.H2II_density[0]  = tiny_number  * my_fields.density[0];
    my_fields.DI_density[0]    = 2.0 * 3.4e-5 * my_fields.density[0];
    my_fields.DII_density[0]   = tiny_number  * my_fields.density[0];
    my_fields.HDI_density[0]   = tiny_number  * my_fields.density[0];
    my_fields.e_density[0]     = tiny_number  * my_fields.density[0];
    
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
    // Atomic
    my_fields.e_density[0]    = rho * *ne_guess;
    
    my_fields.HI_density[0]    = rho * SphP[target].grHI;  //initialized with HYDROGEN_MASSFRAC
    my_fields.HII_density[0]   = rho * SphP[target].grHII;
    my_fields.HM_density[0]    = rho * SphP[target].grHM;
    
    my_fields.HeI_density[0]   = rho * SphP[target].grHeI;
    my_fields.HeII_density[0]  = rho * SphP[target].grHeII;
    my_fields.HeIII_density[0] = rho * SphP[target].grHeIII;
#endif
    
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
    my_fields.H2I_density[0]  = rho * SphP[target].grH2I;
    my_fields.H2II_density[0] = rho * SphP[target].grH2II;
    gamma = SphP[target].Gamma;
#endif
    
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
    my_fields.DI_density[0]   = rho * SphP[target].grDI;
    my_fields.DII_density[0]  = rho * SphP[target].grDII;
    my_fields.HDI_density[0]  = rho * SphP[target].grHDI;
#endif

    // volumetric heating rate (provide in units like [erg s^-1 cm^-3])
    my_fields.volumetric_heating_rate    = malloc(field_size * sizeof(gr_float));
    my_fields.volumetric_heating_rate[0] = 0.0;

    // specific heating rate (provide in units like [egs s^-1 g^-1]
    my_fields.specific_heating_rate      = malloc(field_size * sizeof(gr_float));
    my_fields.specific_heating_rate[0]   = 0.0;

    // radiative transfer ionization / dissociation rate fields (provide in units like*[1/s])
    my_fields.RT_HI_ionization_rate   = malloc(field_size * sizeof(gr_float));
    my_fields.RT_HeI_ionization_rate  = malloc(field_size * sizeof(gr_float));
    my_fields.RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
    my_fields.RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
    // radiative transfer heating rate field (provide in units like [erg s^-1 cm^-3])
    my_fields.RT_heating_rate         = malloc(field_size * sizeof(gr_float));

    my_fields.RT_HI_ionization_rate[0]   = 0.0;
    my_fields.RT_HeI_ionization_rate[0]  = 0.0;
    my_fields.RT_HeII_ionization_rate[0] = 0.0;
    my_fields.RT_H2_dissociation_rate[0] = 0.0;
    my_fields.RT_heating_rate[0]         = 0.0;
    
    switch(mode) {
        case 0:  //solve chemistry & update values
            if(solve_chemistry(&All.GrackleUnits, &my_fields, dt) == 0) {
                fprintf(stderr, "Error in solve_chemistry.\n");
                endrun(ENDRUNVAL);
            }
            
#if (GRACKLE_CHEMISTRY >  0) // non-tabular
            // Assign variables back
            *ne_guess            = my_fields.e_density[0]     / my_fields.density[0];
            
            SphP[target].grHI    = my_fields.HI_density[0]    / my_fields.density[0];
            SphP[target].grHII   = my_fields.HII_density[0]   / my_fields.density[0];
            SphP[target].grHM    = my_fields.HM_density[0]    / my_fields.density[0];
            
            SphP[target].grHeI   = my_fields.HeI_density[0]   / my_fields.density[0];
            SphP[target].grHeII  = my_fields.HeII_density[0]  / my_fields.density[0];
            SphP[target].grHeIII = my_fields.HeIII_density[0] / my_fields.density[0];
#endif
            
#if (GRACKLE_CHEMISTRY >= 2) // Atomic+(H2+H2I+H2II)
            SphP[target].grH2I   = my_fields.H2I_density[0]   / my_fields.density[0];
            SphP[target].grH2II  = my_fields.H2II_density[0]  / my_fields.density[0];
#endif
            
#if (GRACKLE_CHEMISTRY >= 3) // Atomic+(H2+H2I+H2II)+(DI+DII+HD)
            SphP[target].grDI    = my_fields.DI_density[0]    / my_fields.density[0];
            SphP[target].grDII   = my_fields.DII_density[0]   / my_fields.density[0];
            SphP[target].grHDI   = my_fields.HDI_density[0]   / my_fields.density[0];
#endif
            returnval = my_fields.internal_energy[0];
            break;
            
        case 1:  //cooling time
            if(calculate_cooling_time(&All.GrackleUnits, &my_fields,
                                      &cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            if(calculate_temperature(&All.GrackleUnits, &my_fields,
                                     &temperature) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = temperature;
            break;
        case 3:  //calculate pressure
            if(calculate_pressure(&All.GrackleUnits, &my_fields,
                                  &pressure) == 0) {
                fprintf(stderr, "Error in calculate_temperature.\n");
                endrun(ENDRUNVAL);
            }
            returnval = pressure;
            break;
        case 4:  //calculate gamma
            if(calculate_gamma(&All.GrackleUnits, &my_fields,
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

    // Set initial expansion factor (for internal units).
    // Set expansion factor to 1 for non-cosmological simulation.
    double a_value = 1.0;
    if(All.ComovingIntegrationOn) a_value = All.TimeBegin;
    
    // First, set up the units system.
    // These are conversions from code units to cgs.
    // The following parameter must be kept equal to 0 ALWAYS!!!
    All.GrackleUnits.comoving_coordinates = 0; //All.ComovingIntegrationOn; // 1 if cosmological sim, 0 if not
    All.GrackleUnits.density_units        = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
    All.GrackleUnits.length_units         = All.UnitLength_in_cm / All.HubbleParam;
    All.GrackleUnits.time_units           = All.UnitTime_in_s / All.HubbleParam;
    All.GrackleUnits.velocity_units       = All.UnitVelocity_in_cm_per_s;
    All.GrackleUnits.a_units              = 1.0; // units for the expansion factor
    All.GrackleUnits.a_value              = a_value;
    
    // Second, create a chemistry object for parameters and rate data.
    chemistry_data *my_grackle_data;
    my_grackle_data = malloc(sizeof(chemistry_data));
    if (set_default_chemistry_parameters(my_grackle_data) == 0) {
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
    my_grackle_data->three_body_rate        = 0;
    
#ifdef GRACKLE_OPTS
    my_grackle_data->metal_cooling          = All.MetalCooling;    // metal cooling on
#else
    my_grackle_data->metal_cooling          = 0;                   // metal cooling on
#endif
    my_grackle_data->h2_on_dust             = 0;                   // dust cooling/chemistry off
    my_grackle_data->photoelectric_heating            = 0;
    my_grackle_data->photoelectric_heating_rate       = 8.5e-26; 
    
    // Flag to enable an effective CMB temperature floor. This is implemented by subtracting the value of the cooling rate at TCMB from the total cooling rate. Default: 1.
    my_grackle_data->cmb_temperature_floor  = 1;
    // Flag to enable a UV background. If enabled, the cooling table to be used must be specified with the grackle_data_file parameter. Default: 0.
#ifdef GRACKLE_OPTS
    my_grackle_data->UVbackground           = All.UVBackgroundOn;   // UV background on
#else
    my_grackle_data->UVbackground           = 0;                  // UV background on
#endif
    // Flag to enable Compton heating from an X-ray background following Madau & Efstathiou (1999). Default: 0.
    my_grackle_data->Compton_xray_heating   = 1;
    
    
    // Flag to enable H2 collision-induced emission cooling from Ripamonti & Abel (2004). Default: 0.
    my_grackle_data->cie_cooling                      = 0;
    // Flag to enable H2 cooling attenuation from Ripamonti & Abel (2004). Default: 0
    my_grackle_data->h2_optical_depth_approximation   = 0;
    
    // Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field,
    //    in units of 10-21 erg s-1 cm-2 Hz-1 sr-1. Default: 0.
    my_grackle_data->LWbackground_intensity           = 0;
    // Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
    //    (giving a sawtooth pattern), taken from Haiman & Abel, & Rees (2000). Default: 0.
    my_grackle_data->LWbackground_sawtooth_suppression = 0;
    
    
    /* fixed flags:: */
    
    // Flag to activate the grackle machinery:
    my_grackle_data->use_grackle            = 1;                   // grackle on (duh)
    // Path to the data file containing the metal cooling and UV background tables:
    // (expand any ENVVARs in the filename, like $HOME)
    wordexp_t p;
    wordexp(All.GrackleDataFile, &p, 0);
    if(ThisTask == 0)
        printf("expanded GrackleDataFile from %s to %s\n", All.GrackleDataFile, p.we_wordv[0]);
    strcpy(All.GrackleDataFile, p.we_wordv[0]);
    wordfree(&p);

    my_grackle_data->grackle_data_file      = All.GrackleDataFile; // data file
    // Flag to include radiative cooling and actually update the thermal energy during the
    // chemistry solver. If off, the chemistry species will still be updated. The most
    // common reason to set this to off is to iterate the chemistry network to an equilibrium state. Default: 1.
    my_grackle_data->with_radiative_cooling = 1;                   // cooling on
    // The ratio of specific heats for an ideal gas. A direct calculation for the molecular component is used if primordial_chemistry > 1. Default: 5/3.
    my_grackle_data->Gamma                  = GAMMA;              // our eos set in Config.sh
    // Flag to control which primordial chemistry network is used (set by Config file)
#ifndef GRACKLE_CHEMISTRY
    my_grackle_data->primordial_chemistry = 0;                     // fully tabulated cooling
#else
    my_grackle_data->primordial_chemistry = GRACKLE_CHEMISTRY;
#endif
    
    // Finally, initialize the chemistry object.
    if (initialize_chemistry_data(&All.GrackleUnits) == 0) {
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
    species[0]=HYDROGEN_MASSFRAC;//my_grackle_data->HydrogenFractionByMass;
    species[2]=(1-HYDROGEN_MASSFRAC);//my_grackle_data->HydrogenFractionByMass);
    //species[9]=my_grackle_data->DeuteriumToHydrogenRatio*species[0];
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
             SphP[i].Ne     =DMAX(species[5] ,1.0e-20);
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
