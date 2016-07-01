#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "allvars.h"
#include "proto.h"
#include "kernel.h"

/*! \file timestep.c
 *  routines for assigning new timesteps
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO; these 
 * modifications include the addition of various timestep criteria, the WAKEUP 
 * additions, and various changes of units and variable naming conventions throughout. 
 */

static double dt_displacement = 0;

#ifdef SHEARING_BOX
void calc_shearing_box_pos_offset(void)
{
    Shearing_Box_Pos_Offset = Shearing_Box_Vel_Offset * All.Time;
    while(Shearing_Box_Pos_Offset > boxSize_Y) {Shearing_Box_Pos_Offset -= boxSize_Y;}
}
#endif

void set_cosmo_factors_for_current_time(void)
{
    
    /* These are critical factors used throughout for co-moving integrations. Set them here and
       call THESE, rather than trying to come up with the factors throughout, since that makes debugging a nightmare */
    if(All.ComovingIntegrationOn)
    {
        /* All.cf_atime = a = 1/(1+z), the cosmological scale factor */
        All.cf_atime = All.Time;
        /* All.cf_a2inv is just handy */
        All.cf_a2inv = 1 / (All.Time * All.Time);
        /* All.cf_a3inv * Density_code = Density_physical */
        All.cf_a3inv = 1 / (All.Time * All.Time * All.Time);
        /* Pressure_code/Density_code = All.cf_afac1 * Pressure_physical/Density_physical */
        All.cf_afac1 = 1;
        /* All.cf_afac2 * Pressure_code/Density_code * 1/r_code = Pressure_physical/Density_physical * 1/r_physical */
        All.cf_afac2 = 1 / (All.Time * All.cf_afac1);
        /* All.cf_afac3 * sqrt(Pressure_code/Density_code) = sqrt(Pressure_phys/Density_phys) = cs_physical */
        All.cf_afac3 = 1 / sqrt(All.cf_afac1);
        /* time units: proper time dt_phys = 1/hubble_function(a) * dz/(1+z) = dlna / hubble_function(a)
            code time unit in comoving is dlna, so dt_phys = dt_code / All.cf_hubble_a   */
        All.cf_hubble_a = hubble_function(All.Time); /* hubble_function(a) = H(a) = H(z) */
        /* dt_code * v_code/r_code = All.cf_hubble_a2 * dt_phys * v_phys/r_phys */
        All.cf_hubble_a2 = All.Time * All.Time * hubble_function(All.Time);
    }
    else
    {
        All.cf_atime = 1;
        All.cf_a2inv = 1;
        All.cf_a3inv = 1;
        All.cf_afac1 = 1;
        All.cf_afac2 = 1;
        All.cf_afac3 = 1;
        All.cf_hubble_a = 1;
        All.cf_hubble_a2 = 1;
    }
}


/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void find_timesteps(void)
{
    CPU_Step[CPU_MISC] += measure_time();
    
    int i, bin, binold, prev, next;
    integertime ti_step, ti_step_old, ti_min;
    double aphys;
    
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin || dt_displacement == 0)
        find_dt_displacement_constraint(All.cf_hubble_a * All.cf_atime * All.cf_atime);
    
    
#ifdef DIVBCLEANING_DEDNER
    /* need to calculate the global fastest wave speed to manage the damping terms stably */
    if((All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)||(All.FastestWaveSpeed == 0))
    {
        double fastwavespeed = 0.0;
        double fastwavedecay = 0.0;
        double fac_magnetic_pressure = All.cf_afac1 / All.cf_atime;
        for(i=0;i<NumPart;i++)
        {
            if(P[i].Type==0)
            {
                double vsig2 = 0.5 * All.cf_afac3 * fabs(SphP[i].MaxSignalVel); // in v_phys units //
                double vsig1 = All.cf_afac3 * sqrt( Particle_effective_soundspeed_i(i)*Particle_effective_soundspeed_i(i) + fac_magnetic_pressure * (Get_Particle_BField(i,0)*Get_Particle_BField(i,0)+Get_Particle_BField(i,1)*Get_Particle_BField(i,1)+Get_Particle_BField(i,2)*Get_Particle_BField(i,2)) / SphP[i].Density );
                double vsig0 = DMAX(vsig1,vsig2);

                if(vsig0 > fastwavespeed) fastwavespeed = vsig0; // physical unit
                double hsig0 = Get_Particle_Size(i) * All.cf_atime; // physical unit
                if(vsig0/hsig0 > fastwavedecay) fastwavedecay = vsig0 / hsig0; // physical unit
            }
        }
        /* if desired, can just do this by domain; otherwise we use an MPI call over all domains to collect */
        double fastwavespeed_max_glob=fastwavespeed;
        MPI_Allreduce(&fastwavespeed, &fastwavespeed_max_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double fastwavedecay_max_glob=fastwavedecay;
        MPI_Allreduce(&fastwavedecay, &fastwavedecay_max_glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        /* now set the variables */
        All.FastestWaveSpeed = fastwavespeed_max_glob;
        All.FastestWaveDecay = fastwavedecay_max_glob;
    }
#endif

#ifdef FORCE_EQUAL_TIMESTEPS
    for(i = FirstActiveParticle, ti_min = TIMEBASE; i >= 0; i = NextActiveParticle[i])
    {
        ti_step = get_timestep(i, &aphys, 0);
        
        if(ti_step < ti_min)
            ti_min = ti_step;
    }
    
    if(ti_min > (dt_displacement / All.Timebase_interval))
        ti_min = (dt_displacement / All.Timebase_interval);
    
    ti_step = TIMEBASE;
    while(ti_step > ti_min)
        ti_step >>= 1;
    
    integertime ti_min_glob;
    
    MPI_Allreduce(&ti_step, &ti_min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    
    
    /* Now assign new timesteps  */
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
#ifdef FORCE_EQUAL_TIMESTEPS
        ti_step = ti_min_glob;
#else
        ti_step = get_timestep(i, &aphys, 0);
#endif
        
        /* make it a power 2 subdivision */
        ti_min = TIMEBASE;
        while(ti_min > ti_step)
            ti_min >>= 1;
        ti_step = ti_min;
        
        bin = get_timestep_bin(ti_step);
        binold = P[i].TimeBin;
        
        if(bin > binold)		/* timestep wants to increase */
        {
            while(TimeBinActive[bin] == 0 && bin > binold)	/* make sure the new step is synchronized */
                bin--;
            
            ti_step = bin ? (((integertime) 1) << bin) : 0;
        }
        
        if(All.Ti_Current >= TIMEBASE)	/* we here finish the last timestep. */
        {
            ti_step = 0;
            bin = 0;
        }
        
        if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
        {
            terminate("we are beyond the end of the timeline");	/* should not happen */
            ti_step = TIMEBASE - All.Ti_Current;
            ti_min = TIMEBASE;
            while(ti_min > ti_step)
                ti_min >>= 1;
            ti_step = ti_min;
        }
        
        if(bin != binold)
        {
            TimeBinCount[binold]--;
            if(P[i].Type == 0)
            {
                TimeBinCountSph[binold]--;
#ifdef SFR
                TimeBinSfr[binold] -= SphP[i].Sfr;
                TimeBinSfr[bin] += SphP[i].Sfr;
#endif
            }
            
            
            prev = PrevInTimeBin[i];
            next = NextInTimeBin[i];
            
            if(FirstInTimeBin[binold] == i)
                FirstInTimeBin[binold] = next;
            if(LastInTimeBin[binold] == i)
                LastInTimeBin[binold] = prev;
            if(prev >= 0)
                NextInTimeBin[prev] = next;
            if(next >= 0)
                PrevInTimeBin[next] = prev;
            
            if(TimeBinCount[bin] > 0)
            {
                PrevInTimeBin[i] = LastInTimeBin[bin];
                NextInTimeBin[LastInTimeBin[bin]] = i;
                NextInTimeBin[i] = -1;
                LastInTimeBin[bin] = i;
            }
            else
            {
                FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
                PrevInTimeBin[i] = NextInTimeBin[i] = -1;
            }
            TimeBinCount[bin]++;
            if(P[i].Type == 0)
                TimeBinCountSph[bin]++;
            
            P[i].TimeBin = bin;
        }
        
#ifndef WAKEUP
        ti_step_old = binold ? (((integertime) 1) << binold) : 0;
#else
        ti_step_old = P[i].dt_step;
#endif
        
        P[i].Ti_begstep += ti_step_old;
        
#if defined(WAKEUP) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        P[i].dt_step = ti_step;
#endif
    }
    
    
    
    
    
#ifdef WAKEUP
    process_wake_ups();
#endif
    
    CPU_Step[CPU_TIMELINE] += measure_time();
}



/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
integertime get_timestep(int p,		/*!< particle index */
                         double *aphys,	/*!< acceleration (physical units) */
                         int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding aphys */ )
{
    double ax, ay, az, ac;
    double csnd = 0, dt = 0, dt_courant = 0, dt_divv = 0;
    integertime ti_step;
    int k; k=0;
#ifdef CHEMCOOL
    double hubble_param;
    
    if(All.ComovingIntegrationOn)
        hubble_param = All.HubbleParam;
    else
        hubble_param = 1.0;
#endif
    
    
    if(flag == 0)
    {
        ax = All.cf_a2inv * P[p].GravAccel[0];
        ay = All.cf_a2inv * P[p].GravAccel[1];
        az = All.cf_a2inv * P[p].GravAccel[2];
        
        
        
        if(P[p].Type == 0)
        {
            ax += SphP[p].HydroAccel[0];
            ay += SphP[p].HydroAccel[1];
            az += SphP[p].HydroAccel[2];
        }
        
        ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
        *aphys = ac;
    }
    else
        ac = *aphys;
    
    if(ac == 0) ac = 1.0e-30;
    
    
    if(flag > 0)
    {
        /* this is the non-standard mode; use timestep to get the maximum acceleration tolerated */
        dt = flag * All.Timebase_interval;
        dt /= All.cf_hubble_a;	/* convert dloga to physical timestep  */
        
        ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * All.ForceSoftening[P[p].Type] / (dt * dt);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * DMAX(PPP[p].AGS_Hsml,All.ForceSoftening[P[p].Type]) / (dt * dt);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if(P[p].Type==0)
            ac = 2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * DMAX(PPP[p].Hsml,All.ForceSoftening[P[p].Type]) / (dt * dt);
#endif
        *aphys = ac;
        return flag;
    }
    
    dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * All.ForceSoftening[P[p].Type] / ac);
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime  * KERNEL_CORE_SIZE * DMAX(PPP[p].AGS_Hsml,All.ForceSoftening[P[p].Type]) / ac);
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
    if(P[p].Type == 0)
    {
        dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * KERNEL_CORE_SIZE * DMAX(PPP[p].Hsml,All.ForceSoftening[P[p].Type]) / ac);
    }
#endif


    
    
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    /* make sure smoothing length of non-gas particles doesn't change too much in one timestep */
    if(P[p].Type > 0)
    {
        if(PPP[p].AGS_Hsml > 1.01*All.ForceSoftening[P[p].Type])
        {
            double divVel = fabs(P[p].Particle_DivVel);
            double tmp_divVel = 0, tmp_ags_h = PPP[p].AGS_Hsml;
            int k; for(k=0;k<3;k++) {tmp_divVel += P[p].Vel[k]*P[p].Vel[k];}
#ifdef GALSF
            if(((P[p].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[p].Type == 2)||(P[p].Type==3))))&&(P[p].Mass>0)) {tmp_ags_h = DMAX(tmp_ags_h , PPP[p].Hsml);}
#endif
            tmp_divVel = sqrt(tmp_divVel) / tmp_ags_h;
            divVel = All.cf_a2inv * DMIN(divVel , tmp_divVel);
            if(divVel != 0)
            {
                dt_divv = 0.25 / divVel;
                if(dt_divv < dt) {dt = dt_divv;}
            }
        }
    }
#endif

    
    
    if((P[p].Type == 0) && (P[p].Mass > 0))
        {
            csnd = 0.5 * SphP[p].MaxSignalVel * All.cf_afac3;
            double L_particle = Get_Particle_Size(p);
            
            dt_courant = All.CourantFac * (L_particle*All.cf_atime) / csnd;
            if(dt_courant < dt) dt = dt_courant;

            double dt_prefac_diffusion;
            dt_prefac_diffusion = 0.5;
#if defined(GALSF) || defined(DIFFUSION_OPTIMIZERS)
            dt_prefac_diffusion = 1.8;
#endif
            
            
            


#ifdef MHD_NON_IDEAL
            {
                double b_grad = 0, b_mag = 0;
                int k; for(k=0;k<3;k++)
                {
                    int k2;
                    for(k2=0;k2<3;k2++)
                    {
                        double tmp_grad = SphP[p].Gradients.B[k][k2];
                        b_grad += tmp_grad * tmp_grad;
                    }
                    double tmp_grad = Get_Particle_BField(p,k);
                    b_mag += tmp_grad * tmp_grad;
                }
                double L_cond_inv = sqrt(b_grad / (1.e-37 + b_mag));
                double L_cond = DMAX(L_particle , 1./(L_cond_inv + 1./L_particle)) * All.cf_atime;
                double diff_coeff = SphP[p].Eta_MHD_OhmicResistivity_Coeff + SphP[p].Eta_MHD_HallEffect_Coeff + SphP[p].Eta_MHD_AmbiPolarDiffusion_Coeff;
                double dt_conduction = dt_prefac_diffusion * L_cond*L_cond / (1.0e-37 + diff_coeff);
                if(dt_conduction < dt) dt = dt_conduction; // normal explicit time-step
            }
#endif
            
            

            
    
    
            
            

#ifdef TURB_DIFFUSION
            {
            }
#endif
            
            
#if defined(DIVBCLEANING_DEDNER) 
            double fac_magnetic_pressure = All.cf_afac1 / All.cf_atime;
            double phi_b_units = Get_Particle_PhiField(p) / (All.cf_afac3 * All.cf_atime * SphP[p].MaxSignalVel);
            double vsig1 = All.cf_afac3 * sqrt( Particle_effective_soundspeed_i(p)*Particle_effective_soundspeed_i(p) +
                    fac_magnetic_pressure * (Get_Particle_BField(p,0)*Get_Particle_BField(p,0) +
                                             Get_Particle_BField(p,1)*Get_Particle_BField(p,1)+
                                             Get_Particle_BField(p,2)*Get_Particle_BField(p,2) +
                                             phi_b_units*phi_b_units) / SphP[p].Density );

            dt_courant = 0.8 * All.CourantFac * (All.cf_atime*L_particle) / vsig1; // 2.0 factor may be added (PFH) //
            if(dt_courant < dt) {dt = dt_courant;}
#endif
            
            /* make sure that the velocity divergence does not imply a too large change of density or kernel length in the step */
            double divVel = P[p].Particle_DivVel;
            if(divVel != 0)
            {
                dt_divv = 1.5 / fabs(All.cf_a2inv * divVel);
                if(dt_divv < dt) {dt = dt_divv;}
            }
            
            
            
            
        } // closes if(P[p].Type == 0) [gas particle check] //
    
    
    
    
    // add a 'stellar evolution timescale' criterion to the timestep, to prevent too-large jumps in feedback //
#if defined(YOUNGSTARWINDDRIVING) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    if(((P[p].Type == 4)||((All.ComovingIntegrationOn==0)&&((P[p].Type == 2)||(P[p].Type==3))))&&(P[p].Mass>0))
    {
        double star_age = evaluate_stellar_age_Gyr(P[p].StellarAge);
        double dt_stellar_evol;
        if(star_age<0.1)
        {
            dt_stellar_evol = DMAX(2.0e-4, star_age/250.); // restrict to small steps for young stars //
        } else {
            dt_stellar_evol = star_age/10.;
        }
        // PFH: temporarily modifying the terms above while Marcel studies them: turns out not to be necessary to use as strict a mass-dependent timestep, so faster to comment out //
        //dt_stellar_evol /= ( 1. + 0.1*(P[p].Mass*All.UnitMass_in_g)/(1.0e4*1.989e33) ); // multiply by inverse particle mass, since goal is to prevent too much energy in one time //
        double mcorr = 0.1 * (P[p].Mass*All.UnitMass_in_g/All.HubbleParam)/(1.0e4*1.989e33);
        if(mcorr < 1 && mcorr > 0) {dt_stellar_evol /= mcorr;}
        if(dt_stellar_evol < 1.e-6) {dt_stellar_evol = 1.e-6;}
        dt_stellar_evol /= (0.001*All.UnitTime_in_Megayears/All.HubbleParam); // convert to code units //
        if(dt_stellar_evol>0)
            if(dt_stellar_evol<dt)
                dt = dt_stellar_evol;
    }
#endif

#ifdef GALSF_FB_LUPI
    if(P[p].Type == 4 && All.SNeIIFraction == 0)
    {
        double star_age = StellarLifeTime(SNII_MAX_MASS,P[p].Metallicity[0]/LUPI_SOLAR_MET);
	double current_age = evaluate_stellar_age_Gyr(P[p].StellarAge)*1.e3;
        double dt_stellar_evol = star_age / All.UnitTime_in_Megayears * All.HubbleParam * 0.1; //1/10th of the stellar lifetime
	if(dt_stellar_evol>0 && current_age< SNI_MIN_AGE)
            if(dt_stellar_evol<dt)
		dt = dt_stellar_evol;
    }
#endif
    
#if defined(GRACKLE_CHEMISTRY) && (GRACKLE_CHEMISTRY>0)
    if(P[p].Type==0)
#ifdef GALSF_FB_LUPI
      if(SphP[p].DelayTimeCoolingSNe<=0.)
#endif
      {
	double dt_cool = CallGrackle(SphP[p].InternalEnergyPred, SphP[p].Density, 0, &(SphP[p].Ne), p, 1);
       	dt_cool *= All.HubbleParam;
        if(dt_cool>0)
            if(dt_cool<0.1*dt)
                dt=dt_cool;
      }
#endif
    
    

    
    /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected, All.cf_hubble_a=1. */
    dt *= All.cf_hubble_a;
    
#ifdef ONLY_PM
    dt = All.MaxSizeTimestep;
#endif
    
    
    
    if(dt >= All.MaxSizeTimestep)
        dt = All.MaxSizeTimestep;
    
    
    if(dt >= dt_displacement)
        dt = dt_displacement;
    
    
    if(dt < All.MinSizeTimestep)
    {
#ifdef STOP_WHEN_BELOW_MINTIMESTEP
        printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
        
        if(P[p].Type == 0)
        {
#ifndef LONGIDS
            printf
            ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
             (int) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
             PPP[p].Hsml, csnd,
             sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) *
             All.cf_hubble_a, All.SofteningTable[P[p].Type]);
#else
            printf
            ("Part-ID=%llu  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
             (MyIDType) P[p].ID, dt, dt_courant * All.cf_hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2],
             PPP[p].Hsml, csnd,
             sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.SofteningTable[P[p].Type] / ac) *
             All.cf_hubble_a, All.SofteningTable[P[p].Type]);
#endif // ndef LONGIDS
            
            
        }
        else // if(P[p].Type == 0)
        {
#ifndef LONGIDS
            printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
                   P[p].Pos[2]);
#else
            printf("Part-ID=%llu  dt=%g ac=%g xyz=(%g|%g|%g)\n", (MyIDType) P[p].ID, dt, ac, P[p].Pos[0],
                   P[p].Pos[1], P[p].Pos[2]);
#endif // ndef LONGIDS
        }
        fflush(stdout);
        fprintf(stderr, "\n @ fflush \n");
        endrun(888);
#endif // STOP_WHEN_BELOW_MINTIMESTEP
        dt = All.MinSizeTimestep;
    }
    
    ti_step = (integertime) (dt / All.Timebase_interval);
#ifndef STOP_WHEN_BELOW_MINTIMESTEP
    if(ti_step<=1) ti_step=2;
#endif
    
    
    if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
        printf("\nError: A timestep of size zero was assigned on the integer timeline, no here!!!\n"
               "We better stop.\n"
               "Task=%d Part-ID=%llu dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
               ThisTask, (unsigned long long) P[p].ID, dt, dt_courant, dt_divv, dt_displacement,
               All.Timebase_interval, ti_step, ac,
               P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].GravAccel[0], P[p].GravAccel[1],
               P[p].GravAccel[2]);
        
        fflush(stdout);
        endrun(818);
    }
    
    return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
    int i, type;
    int count[6];
    long long count_sum[6];
    double v[6], v_sum[6], mim[6], min_mass[6];
    double dt, dmean, asmth = 0;
    
    dt_displacement = All.MaxSizeTimestep;
    
    if(All.ComovingIntegrationOn)
    {
        for(type = 0; type < 6; type++)
        {
            count[type] = 0;
            v[type] = 0;
            mim[type] = 1.0e30;
        }
        
        for(i = 0; i < NumPart; i++)
        {
            v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
            if(P[i].Mass > 0)
            {
                if(mim[P[i].Type] > P[i].Mass)
                    mim[P[i].Type] = P[i].Mass;
            }
            count[P[i].Type]++;
        }
        
        MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        sumup_large_ints(6, count, count_sum);
        
#ifdef SFR
        /* add star and gas particles together to treat them on equal footing, using the original gas particle spacing. */
        v_sum[0] += v_sum[4];
        count_sum[0] += count_sum[4];
        v_sum[4] = v_sum[0];
        count_sum[4] = count_sum[0];
#endif
        
        for(type = 0; type < 6; type++)
        {
            if(count_sum[type] > 0)
            {
#ifdef GALSF
                if(type == 0 || type == 4)
#else
                if(type == 0)
#endif
                    dmean = pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);
                else
                    dmean = pow(min_mass[type] / ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G)), 1.0 / 3);
                
                dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);
                
                
                if(ThisTask == 0)
                    printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
                           type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);
                
                if(dt < dt_displacement)
                    dt_displacement = dt;
            }
        }
        
        if(ThisTask == 0)
            printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}



int get_timestep_bin(integertime ti_step)
{
    int bin = -1;
    
    if(ti_step == 0)
        return 0;
    
    if(ti_step == 1)
        terminate("time-step of integer size 1 not allowed\n");
    
    while(ti_step)
    {
        bin++;
        ti_step >>= 1;
    }
    
    return bin;
}





#ifdef WAKEUP
void process_wake_ups(void)
{
    int i, n, dt_bin;
    int ti_next_for_bin, ti_next_kick, ti_next_kick_global, max_time_bin_active;
    int bin, binold, prev, next;
    long long ntot;
    
    /* find the next kick time */
    for(n = 0, ti_next_kick = TIMEBASE; n < TIMEBINS; n++)
    {
        if(TimeBinCount[n])
        {
            if(n > 0)
            {
                dt_bin = (((integertime) 1) << n);
                ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin;	/* next kick time for this timebin */
            }
            else
            {
                dt_bin = 0;
                ti_next_for_bin = All.Ti_Current;
            }
            
            if(ti_next_for_bin < ti_next_kick)
                ti_next_kick = ti_next_for_bin;
        }
    }
    
    MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    
    if(ThisTask == 0)
        printf("predicting next timestep: %g\n", (ti_next_kick_global - All.Ti_Current) * All.Timebase_interval);
    
    max_time_bin_active = 0;
    /* get the highest bin, that is active next time */
    for(n = 0; n < TIMEBINS; n++)
    {
        dt_bin = (((integertime) 1) << n);
        
        if((ti_next_kick_global % dt_bin) == 0)
            max_time_bin_active = n;
    }
    
    /* move the particle on the highest bin, that is active in the next timestep and that is lower than its last timebin */
    bin = 0;
    for(n = 0; n < TIMEBINS; n++)
    {
        if(TimeBinCount[n] > 0)
        {
            bin = n;
            break;
        }
    }
    n = 0;
    
    for(i = 0; i < NumPart; i++)
    {
        if(P[i].Type != 0)
            continue;
        
        if(P[i].Mass <= 0)
            continue;
        
        if(!SphP[i].wakeup)
            continue;
        
        binold = P[i].TimeBin;
        if(TimeBinActive[binold])
            continue;
        
        bin = max_time_bin_active < binold ? max_time_bin_active : binold;
        
        if(bin != binold)
        {
            integertime dt_0 = P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0;
            integertime tstart = P[i].Ti_begstep + dt_0;
            integertime t_2 = P[i].Ti_current;
            if(t_2 > tstart) {tstart = t_2;}
            integertime tend = All.Ti_Current;

            TimeBinCount[binold]--;
            if(P[i].Type == 0)
                TimeBinCountSph[binold]--;
            
            prev = PrevInTimeBin[i];
            next = NextInTimeBin[i];
            
            if(FirstInTimeBin[binold] == i)
                FirstInTimeBin[binold] = next;
            if(LastInTimeBin[binold] == i)
                LastInTimeBin[binold] = prev;
            if(prev >= 0)
                NextInTimeBin[prev] = next;
            if(next >= 0)
                PrevInTimeBin[next] = prev;
            
            if(TimeBinCount[bin] > 0)
            {
                PrevInTimeBin[i] = LastInTimeBin[bin];
                NextInTimeBin[LastInTimeBin[bin]] = i;
                NextInTimeBin[i] = -1;
                LastInTimeBin[bin] = i;
            }
            else
            {
                FirstInTimeBin[bin] = LastInTimeBin[bin] = i;
                PrevInTimeBin[i] = NextInTimeBin[i] = -1;
            }
            TimeBinCount[bin]++;
            if(P[i].Type == 0)
                TimeBinCountSph[bin]++;
            
            P[i].TimeBin = bin;
            
            if(TimeBinActive[bin])
                NumForceUpdate++;
                        
            n++;

            /* reverse part of the last second-half kick this particle received 
                (to correct it back to its new active time) */
            if(tend < tstart)
            {
                do_the_kick(i, tstart, tend, P[i].Ti_current, 1);
                set_predicted_sph_quantities_for_extra_physics(i);
            }
            P[i].Ti_begstep = All.Ti_Current;
            P[i].dt_step = bin ? (((integertime) 1) << bin) : 0;
            if(P[i].Ti_current < All.Ti_Current) {P[i].Ti_current=All.Ti_Current;}
            
        }
    }
    
    sumup_large_ints(1, &n, &ntot);
    if(ThisTask == 0)
    {
        if(ntot > 0) {printf("%d%09d particles woken up.\n", (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
    }
}
#endif
