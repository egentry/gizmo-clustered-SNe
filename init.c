#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"

/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly initializing 
 * new/modified variables, as needed)
 */

/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the initial gas kernel lengths are determined.
 */
void init(void)
{
    int i, j;
    double a3, atime;
    
#ifdef MAGNETIC
    double a2_fac;
    double gauss2gizmo = All.UnitMagneticField_in_gauss / sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam);
    /* NOTE: we will always work -internally- in code units where MU_0 = 1; hence the 4pi here;
        [much simpler, but be sure of your conversions!] */
#endif
    
    
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();    
    
    if(RestartFlag == 3 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf("Need to give the snapshot number if FLAG_NOT_IN_PUBLIC_CODE/FLAG_NOT_IN_PUBLIC_CODE is selected for output\n");
        endrun(0);
    }
    
    if(RestartFlag == 4 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf("Need to give the snapshot number if snapshot should be converted\n");
        endrun(0);
    }
    
    if(RestartFlag == 5 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf
            ("Need to give the snapshot number if power spectrum and two-point correlation function should be calculated\n");
        endrun(0);
    }
    
    if(RestartFlag == 6 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf
            ("Need to give the snapshot number if velocity power spectrum for the gas cells should be calculated\n");
        endrun(0);
    }
    
#ifdef GRACKLE_FIX_TEMPERATURE
    if(RestartFlag == 7 && RestartSnapNum < 0)
    {
        if(ThisTask == 0)
            printf
            ("Need to give the snapshot number if temperature for gas should be calculated\n");
        endrun(0);
    }
#endif
    
    switch (All.ICFormat)
    {
        case 1:
        case 2:
        case 3:
        case 4:
            if(RestartFlag >= 2 && RestartSnapNum >= 0)
            {
                char fname[1000];
                
                if(All.NumFilesPerSnapshot > 1)
                    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
                            RestartSnapNum);
                else
                    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
                
                read_ic(fname);
                
            }
            else
            {
                read_ic(All.InitCondFile);
            }
            break;
            
        default:
            if(ThisTask == 0)
                printf("ICFormat=%d not supported.\n", All.ICFormat);
            endrun(0);
    }
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();
    
    
    
#if defined(COOLING)
    IonizeParams();
#endif
    
    if(All.ComovingIntegrationOn)
    {
        All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
        All.Ti_Current = 0;
        a3 = All.Time * All.Time * All.Time;
        atime = All.Time;
#ifdef MAGNETIC
        a2_fac = (All.Time * All.Time);
#endif
    }
    else
    {
        All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
        All.Ti_Current = 0;
        a3 = 1;
        atime = 1;
#ifdef MAGNETIC
        a2_fac = 1;
#endif
    }
        
    set_softenings();
    
    All.NumCurrentTiStep = 0;	/* setup some counters */
    All.SnapshotFileCount = 0;
    if(RestartFlag == 2)
    {
        if(RestartSnapNum < 0)
        {
            char *underscore = strrchr(All.InitCondFile, '_');
            if(!underscore)
            {
                char buf[1000];
                sprintf(buf, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n",
                        All.InitCondFile);
                terminate(buf);
            }
            else
                All.SnapshotFileCount = atoi(underscore + 1) + 1;
        }
        else
            All.SnapshotFileCount = RestartSnapNum + 1;
    }
    
    
    All.TotNumOfForces = 0;
    All.TopNodeAllocFactor = 0.008; /* this will start from a low value and be iteratively increased until it is well-behaved */
    All.TreeAllocFactor = 0.45; /* this will also iteratively increase to fit the particle distribution */
    /* To construct the BH-tree for N particles, somewhat less than N
     internal tree-nodes are necessary for ‘normal’ particle distributions. 
     TreeAllocFactor sets the number of internal tree-nodes allocated in units of the particle number. 
     By experience, space for ≃ 0.65N internal nodes is usually fully sufficient for typical clustered 
     particle distributions, so a value of 0.7 should put you on the safe side. If the employed particle 
     number per processor is very small (less than a thousand or so), or if there are many particle pairs 
     with identical or nearly identical coordinates, a higher value may be required. Since the number of 
     particles on a given processor may be higher by a factor PartAllocFactor than the average particle 
     number, the total amount of memory requested for the BH tree on a single processor scales proportional 
     to PartAllocFactor*TreeAllocFactor. */
    
    
#ifdef GRACKLE_FIX_TEMPERATURE
    if(RestartFlag == 7)
    {
        compute_temperature();
        sprintf(All.SnapshotFileBase, "%s_temp", All.SnapshotFileBase);
        printf("RestartSnapNum %d\n", RestartSnapNum);
        savepositions(RestartSnapNum);
        endrun(0);
    }
#endif
    
#ifdef PERIODIC
    if(All.ComovingIntegrationOn) check_omega();
#endif
    
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif
    
    for(i = 0; i < GRAVCOSTLEVELS; i++)
        All.LevelToTimeBin[i] = 0;
    
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < GRAVCOSTLEVELS; j++)
            P[i].GravCost[j] = 0;
    
    
    
    if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
        for(i = 0; i < NumPart; i++)
            for(j = 0; j < 3; j++)
	    {
                P[i].Vel[j] *= sqrt(All.Time) * All.Time;
	    }
    }
    
    
    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
        for(j = 0; j < 3; j++)
            P[i].GravAccel[j] = 0;
        
        /* DISTORTION PARTICLE SETUP */
        
#ifdef KEEP_DM_HSML_AS_GUESS
        if(RestartFlag != 1)
            P[i].DM_Hsml = -1;
#endif
        
        P[i].Ti_begstep = 0;
        P[i].Ti_current = 0;
        P[i].TimeBin = 0;
        
        if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
            P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */
        
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY)    
        P[i].Potential = 0;
#endif
#ifdef GALSF
        if(RestartFlag == 0)
        {
            P[i].StellarAge = 0;
        }
#endif
        
        if(RestartFlag != 1)
        {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES)
            P[i].DensAroundStar = 0;
            P[i].GradRho[0]=0;
            P[i].GradRho[1]=0;
            P[i].GradRho[2]=1;
#endif
#ifdef GALSF_FB_LUPI
	    P[i].SNe_ThisTimeStep = 0;
	    P[i].MassYield_ThisTimeStep = 0;
	    P[i].MassLoss_ThisTimeStep = 0;
	    P[i].MetalYield_ThisTimeStep[0] = 0; //Z is always evolved
	    if(NUM_METAL_SPECIES>1)
	    {
		P[i].MetalYield_ThisTimeStep[1] = 0; //Oxygen
		P[i].MetalYield_ThisTimeStep[2] = 0; //Iron
	    }
	    P[i].WeightNorm[0] = 0;
	    P[i].WeightNorm[1] = 0;
#endif

        }
        
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        if(RestartFlag == 0)
        {
            P[i].StellarAge = -2.0 * All.InitStellarAgeinGyr / (All.UnitTime_in_Megayears*0.001) * get_random_number(P[i].ID + 3);
        }
#endif

#ifdef GALSF_FB_LUPI
	if(RestartFlag == 0)
	{
	    P[i].StellarAge = MAX_REAL_NUMBER;
	    P[i].StellarInitMass = 0;
	}
#endif

#ifdef BH_LUPI
	if(RestartFlag == 0)
	    P[i].BH_StoredEnergy = 0;

	if(RestartFlag != 1)
	{
	    P[i].BH_Luminosity = P[i].BH_AccretionRate = 0;
	    P[i].BH_DeltaPos[0] = 0;
	    P[i].BH_DeltaPos[1] = 0;
	    P[i].BH_DeltaPos[2] = 0;
	    P[i].BH_DeltaVel[0] = 0;
	    P[i].BH_DeltaVel[1] = 0;
	    P[i].BH_DeltaVel[2] = 0;
	    P[i].BH_GasVelocity[0] = 0;
	    P[i].BH_GasVelocity[1] = 0;
	    P[i].BH_GasVelocity[2] = 0;
	}
#endif

        
        
        
        
        
        
#ifdef GRACKLE_OPTS
            if(RestartFlag == 0 && P[i].Type!=1 && P[i].Type!=5)
		P[i].Metallicity[0]=All.InitMetallicityinSolar*GENTRY_SOLAR_MET;
#endif
        
        
    }
    
    
    for(i = 0; i < TIMEBINS; i++)
        TimeBinActive[i] = 1;
    
    reconstruct_timebins();
    
        
    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        
        for(j = 0; j < 3; j++)
        {
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
            //SphP[i].dMomentum[j] = 0;//manifest-indiv-timestep-debug//
        }
        
        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        P[i].Particle_DivVel = 0;
        SphP[i].ConditionNumber = 1;
        SphP[i].DtInternalEnergy = 0;
#ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        SphP[i].MaxKineticEnergyNgb = 0;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif
        
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[i].AGS_zeta = 0;
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        PPP[i].AGS_Hsml = PPP[i].Hsml;
#endif
#endif
        
#ifdef MHD_NON_IDEAL
        SphP[i].Eta_MHD_OhmicResistivity_Coeff = 0;
        SphP[i].Eta_MHD_HallEffect_Coeff = 0;
        SphP[i].Eta_MHD_AmbiPolarDiffusion_Coeff = 0;
#endif
        
        
#ifdef TURB_DIFFUSION
        SphP[i].TD_DiffCoeff = 0;
#endif
        
        if(RestartFlag == 0)
        {
#ifndef READ_HSML
            PPP[i].Hsml = 0;
#endif
            SphP[i].Density = -1;
#ifdef COOLING
            SphP[i].Ne = 1.0;
#endif
        }
#ifdef GALSF_FB_LUPI
	if(RestartFlag == 0) SphP[i].DelayTimeCoolingSNe = 0;
#endif

#ifdef GALSF
        SphP[i].Sfr = 0;
#endif
#ifdef MAGNETIC
#if defined B_SET_IN_PARAMS
        if(RestartFlag == 0)
        {			/* Set only when starting from ICs */
            SphP[i].B[0]=SphP[i].BPred[0] = All.BiniX;
            SphP[i].B[1]=SphP[i].BPred[1] = All.BiniY;
            SphP[i].B[2]=SphP[i].BPred[2] = All.BiniZ;
        }
#endif /*B_SET_IN_PARAMS*/
        for(j = 0; j < 3; j++)
        {
            SphP[i].BPred[j] *= a2_fac * gauss2gizmo;
            SphP[i].B[j] = SphP[i].BPred[j];
        }
#if defined(SPH_TP12_ARTIFICIAL_RESISTIVITY)
        SphP[i].Balpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
        SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif
#ifdef SPHAV_CD10_FLAG_NOT_IN_PUBLIC_CODE_SWITCH
        SphP[i].alpha = 0.0;
#endif
    }
    
#ifndef SHEARING_BOX
#if (NUMDIMS==2)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[2] = 0;
        //P[i].Vel[2] = 0; // this should be set in the ICs, not here //
        
        P[i].GravAccel[2] = 0;
        
        if(P[i].Type == 0)
        {
            SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[2] = 0;
        }
    }
#endif
#endif
    
#if (NUMDIMS==1)
    for(i = 0; i < NumPart; i++)
    {
        P[i].Pos[1] = P[i].Pos[2] = 0;
        //P[i].Vel[1] = P[i].Vel[2] = 0; // this should be set in the ICs, not here //
        
        P[i].GravAccel[1] = P[i].GravAccel[2] = 0;
        
        if(P[i].Type == 0)
        {
            SphP[i].VelPred[1] = SphP[i].VelPred[2] = 0;
            SphP[i].HydroAccel[1] = SphP[i].HydroAccel[2] = 0;
        }
    }
#endif
    
#ifdef ASSIGN_NEW_IDS
    assign_unique_ids();
#endif
    /* assign other ID parameters needed */
    for(i = 0; i < NumPart; i++) {P[i].ID_child_number = 0; P[i].ID_generation = 0;}
    
#ifdef TEST_FOR_IDUNIQUENESS
    test_id_uniqueness();
#endif
    
    Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */
    
    TreeReconstructFlag = 1;
    
    
#ifdef SHIFT_BY_HALF_BOX
    for(i = 0; i < NumPart; i++)
        for(j = 0; j < 3; j++)
        {
            double boxtmp = 0;
            if(j==0) {boxtmp = boxSize_X;}
            if(j==1) {boxtmp = boxSize_Y;}
            if(j==2) {boxtmp = boxSize_Z;}
            P[i].Pos[j] += 0.5 * boxtmp;
        }
#endif
    
    
    Gas_split = 0;
#ifdef GALSF
    Stars_converted = 0;
#endif
    domain_Decomposition(0, 0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */
    
    set_softenings();
    
    /* will build tree */
    ngb_treebuild();
    
    All.Ti_Current = 0;
    
    if(RestartFlag != 3 && RestartFlag != 5)
        setup_smoothinglengths();
    
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    if(RestartFlag != 3 && RestartFlag != 5)
        ags_setup_smoothinglengths();
#endif
    
    

    
    /* HELLO! This here is where you should insert custom code for hard-wiring the ICs of various test problems */

    
    
    density();

#if (GRACKLE_CHEMISTRY>=1) || defined(KROME)
    if(RestartFlag == 0) init_species();
#endif


    for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        
        // re-match the predicted and initial velocities and B-field values, just to be sure //
        for(j=0;j<3;j++) SphP[i].VelPred[j]=P[i].Vel[j];
#ifdef MAGNETIC
        for(j=0;j<3;j++) {SphP[i].B[j] = SphP[i].BPred[j] * P[i].Mass / SphP[i].Density;} // convert to the conserved unit V*B //
        for(j=0;j<3;j++) {SphP[i].BPred[j]=SphP[i].B[j]; SphP[i].DtB[j]=0;}
#endif
        //SphP[i].dInternalEnergy = 0;//manifest-indiv-timestep-debug//
        SphP[i].DtInternalEnergy = 0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[i].dMass = 0;
        SphP[i].DtMass = 0;
        SphP[i].MassTrue = P[i].Mass;
        for(j=0;j<3;j++) SphP[i].GravWorkTerm[j] = 0;
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        PPPZ[i].AGS_zeta = 0;
#endif
        
#ifdef GRACKLE
        if(RestartFlag == 0)
        {
#if (GRACKLE_CHEMISTRY >= 1)
            SphP[i].grHI    = HYDROGEN_MASSFRAC;
            SphP[i].grHII   = 1.0e-20;
            SphP[i].grHM    = 1.0e-20;
            SphP[i].grHeI   = 1.0 - HYDROGEN_MASSFRAC;
            SphP[i].grHeII  = 1.0e-20;
            SphP[i].grHeIII = 1.0e-20;
#endif
#if (GRACKLE_CHEMISTRY >= 2)
            SphP[i].grH2I   = 1.0e-20;
            SphP[i].grH2II  = 1.0e-20;
#endif
#if (GRACKLE_CHEMISTRY >= 3)
            SphP[i].grDI    = 2.0 * 3.4e-5;
            SphP[i].grDII   = 1.0e-20;
            SphP[i].grHDI   = 1.0e-20;
#endif
        }
#endif
        
    }
    
    
    /* we should define the maximum and minimum particle masses 
        below/above which particles are merged/split */
    if(RestartFlag != 1)
    {
        double mass_min = MAX_REAL_NUMBER;
        double mass_max = -MAX_REAL_NUMBER;
        for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
        {
            if(P[i].Mass > mass_max) mass_max = P[i].Mass;
            if(P[i].Mass < mass_min) mass_min = P[i].Mass;
        }
        /* broadcast this and get the min and max values over all processors */
        double mpi_mass_min,mpi_mass_max;
        MPI_Allreduce(&mass_min, &mpi_mass_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_max, &mpi_mass_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        All.MinMassForParticleMerger = 0.49 * mpi_mass_min;
#ifdef GALSF_GENERATIONS
        All.MinMassForParticleMerger /= (float)GALSF_GENERATIONS;
#endif
        /* All.MaxMassForParticleSplit  = 5.01 * mpi_mass_max; */
        All.MaxMassForParticleSplit  = 3.01 * mpi_mass_max;
    }
    
    
    
    
    if(RestartFlag == 3)
    {
        
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(ThisTask == 0)
            printf("*ADAPTIVE_GRAVSOFT_FORALL* Computation of softening lengths... \n");
        ags_setup_smoothinglengths();
        if(ThisTask == 0)
            printf("*ADAPTIVE_GRAVSOFT_FORALL* Computation of softening lengths done. \n");
#endif
        
        endrun(0);
    }
    
    
    
    if(RestartFlag == 6)
    {
        endrun(0);
    }
    
    
    if(RestartFlag == 4)
    {
        All.Time = All.TimeBegin = header.time;
        sprintf(All.SnapshotFileBase, "%s_converted", All.SnapshotFileBase);
        if(ThisTask == 0)
            printf("Start writing file %s\n", All.SnapshotFileBase);
        printf("RestartSnapNum %d\n", RestartSnapNum);
        
        All.TopNodeAllocFactor = 0.008;
        
        savepositions(RestartSnapNum);
        endrun(0);
    }
}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
#ifdef PERIODIC
void check_omega(void)
{
    double mass = 0, masstot, omega;
    int i;
    
    for(i = 0; i < NumPart; i++)
        mass += P[i].Mass;
    
    MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    omega = masstot / (boxSize_X*boxSize_Y*boxSize_Z) / (3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G));
    
    //if(fabs(omega - All.Omega0) > 1.0e-3)
    // because of how we set up these ICs, allow a little more generous tolerance
    if(fabs(omega - All.Omega0) > 1.0e-2)
    {
        if(ThisTask == 0)
        {
            printf("\n\nI've found something odd!\n");
            printf
            ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
             omega, All.Omega0);
            printf("\nI better stop.\n");
            
            fflush(stdout);
        }
        endrun(1);
    }
}
#endif


/*! This function is used to find an initial kernel length (what used to be called the 
 *  'smoothing length' for SPH, but is just the kernel size for the mesh-free methods) for each gas
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the kernel length is provided to the function density(), which will
 *  then iterate if needed to find the right kernel length.
 */
void setup_smoothinglengths(void)
{
    int i, no, p;
    if((RestartFlag == 0)||(RestartFlag==2)) // best for stability if we re-calc Hsml for snapshot restarts //
    {
#if defined(DO_DENSITY_AROUND_STAR_PARTICLES) || defined(FLAG_NOT_IN_PUBLIC_CODE)
        for(i = 0; i < NumPart; i++)
#else
            for(i = 0; i < N_gas; i++)
#endif
            {
                no = Father[i];
                
                while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    
                    if(p < 0)
                        break;
                    
                    no = p;
                }
                
                if((RestartFlag == 0)||(P[i].Type != 0)) // if Restartflag==2, use the saved Hsml of the gas as initial guess //
                {
                    
#ifndef READ_HSML
#if NUMDIMS == 3
                    PPP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.333333) * Nodes[no].len;
#endif
#if NUMDIMS == 2
                    PPP[i].Hsml = pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 0.5) * Nodes[no].len;
#endif
#if NUMDIMS == 1
                    PPP[i].Hsml = All.DesNumNgb * (P[i].Mass / Nodes[no].u.d.mass) * Nodes[no].len;
#endif
#ifndef NOGRAVITY
                    if(All.SofteningTable[0] != 0)
                    {
                        if((PPP[i].Hsml>100.*All.SofteningTable[0])||(PPP[i].Hsml<=0.01*All.SofteningTable[0])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                            PPP[i].Hsml = All.SofteningTable[0];
                    }
#else
                    if((Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0)) PPP[i].Hsml = 1.0;
#endif
#endif // READ_HSML
                } // closes if((RestartFlag == 0)||(P[i].Type != 0))
            }
    }
    
    
    
 
    density();    
}


void assign_unique_ids(void)
{
    int i, *numpartlist;
    MyIDType idfirst;
    
    numpartlist = (int *) mymalloc("numpartlist", NTask * sizeof(int));
    
    MPI_Allgather(&NumPart, 1, MPI_INT, numpartlist, 1, MPI_INT, MPI_COMM_WORLD);
    
    idfirst = 1;
    
    for(i = 0; i < ThisTask; i++)
        idfirst += numpartlist[i];
    
    for(i = 0; i < NumPart; i++)
    {
        P[i].ID = idfirst;
        idfirst++;
    }
    
    myfree(numpartlist);
}


#ifdef ADAPTIVE_GRAVSOFT_FORALL
void ags_setup_smoothinglengths(void)
{
    int i, no, p;
    if(RestartFlag == 0 || RestartFlag == 2)
    {
        for(i = 0; i < NumPart; i++)
        {
            P[i].Particle_DivVel = 0;
            PPPZ[i].AGS_zeta = 0;
            if(P[i].Type > 0)
            {
                no = Father[i];
                while(10 * All.AGS_DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
                {
                    p = Nodes[no].u.d.father;
                    if(p < 0)
                        break;
                    no = p;
                }
                PPP[i].AGS_Hsml = pow(1.0/NORM_COEFF * All.AGS_DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0/NUMDIMS) * Nodes[no].len;
                if(All.SofteningTable[P[i].Type] != 0)
                {
                    if((PPP[i].AGS_Hsml>1000.*All.SofteningTable[P[i].Type])||(PPP[i].AGS_Hsml<=0.01*All.SofteningTable[P[i].Type])||(Nodes[no].u.d.mass<=0)||(Nodes[no].len<=0))
                        PPP[i].AGS_Hsml = All.SofteningTable[P[i].Type];
                }
            } else {
                PPP[i].AGS_Hsml = PPP[i].Hsml;
            }
        }
    }
    ags_density();
}
#endif // ADAPTIVE_GRAVSOFT_FORALL




void test_id_uniqueness(void)
{
    double t0, t1;
#ifndef BND_PARTICLES
    int i;
    MyIDType *ids, *ids_first;
#endif
    
    if(ThisTask == 0)
    {
        printf("Testing ID uniqueness...\n");
        fflush(stdout);
    }
    
    if(NumPart == 0)
    {
        printf("need at least one particle per cpu\n");
        endrun(8);
    }
    
    t0 = my_second();
    
#ifndef BND_PARTICLES
    ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
    ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));
    
    for(i = 0; i < NumPart; i++)
        ids[i] = P[i].ID;
    
#ifdef ALTERNATIVE_PSORT
    init_sort_ID(ids, NumPart);
#else
    parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);
#endif
    
    for(i = 1; i < NumPart; i++)
        if(ids[i] == ids[i - 1])
        {
#ifdef LONGIDS
            printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
                   (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);
            
#else
            printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
            endrun(12);
        }
    
    MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);
    
    if(ThisTask < NTask - 1)
        if(ids[NumPart - 1] == ids_first[ThisTask + 1])
        {
            printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
            endrun(13);
        }
    
    myfree(ids_first);
    myfree(ids);
#endif
    
    t1 = my_second();
    
    if(ThisTask == 0)
    {
        printf("success.  took=%g sec\n", timediff(t0, t1));
        fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
    if(*((MyIDType *) a) < *((MyIDType *) b))
        return -1;
    
    if(*((MyIDType *) a) > *((MyIDType *) b))
        return +1;
    
    return 0;
}

#ifdef GRACKLE_FIX_TEMPERATURE
void compute_temperature()
{
    int i;
    for(i=0;i<N_gas;i++)
    {
       SphP[i].InternalEnergyPred = CallGrackle(SphP[i].InternalEnergy, SphP[i].Density, 0, &(SphP[i].Ne), i, 2);
    }
}
#endif
