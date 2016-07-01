#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>


#include "allvars.h"
#include "proto.h"

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile (except for changing the name of the IC file to the snapshot, 
 * and ensuring the format tag matches it).  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialzed to this
 * value assuming a mean colecular weight either corresponding to complete
 * neutrality, or full ionization.
 */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part (adding/removing read items and changing variable units as necessary)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
static unsigned long FileNr;
static long long *NumPartPerFile;
#endif


void read_ic(char *fname)
{
    int i, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
    double u_init, molecular_weight;
    char buf[500];
    
    CPU_Step[CPU_MISC] += measure_time();
    
#ifdef RESCALEVINI
    if(ThisTask == 0 && RestartFlag == 0)
    {
        fprintf(stdout, "\nRescaling v_ini !\n\n");
        fflush(stdout);
    }
#endif
    
    NumPart = 0;
    N_gas = 0;
    All.TotNumPart = 0;
    
    num_files = find_files(fname);
    
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
    NumPartPerFile = (long long *) mymalloc("NumPartPerFile", num_files * sizeof(long long));
    
    if(ThisTask == 0)
        get_particle_numbers(fname, num_files);
    
    MPI_Bcast(NumPartPerFile, num_files * sizeof(long long), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    
    rest_files = num_files;
    
    while(rest_files > NTask)
    {
        sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
        if(All.ICFormat == 3)
            sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
        FileNr = ThisTask + (rest_files - NTask);
#endif
        
        ngroups = NTask / All.NumFilesWrittenInParallel;
        if((NTask % All.NumFilesWrittenInParallel))
            ngroups++;
        groupMaster = (ThisTask / ngroups) * ngroups;
        
        for(gr = 0; gr < ngroups; gr++)
        {
            if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */
                read_file(buf, ThisTask, ThisTask);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        rest_files -= NTask;
    }
    
    
    if(rest_files > 0)
    {
        distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);
        
        if(num_files > 1)
        {
            sprintf(buf, "%s.%d", fname, filenr);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.%d.hdf5", fname, filenr);
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
            FileNr = filenr;
#endif
        }
        else
        {
            sprintf(buf, "%s", fname);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.hdf5", fname);
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
            FileNr = 0;
#endif
        }
        
        ngroups = rest_files / All.NumFilesWrittenInParallel;
        if((rest_files % All.NumFilesWrittenInParallel))
            ngroups++;
        
        for(gr = 0; gr < ngroups; gr++)
        {
            if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
                read_file(buf, masterTask, lastTask);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    
    
    
    myfree(CommBuffer);
    
    
    if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
    {
        /* this makes sure that masses are initialized in the case that the mass-block
         is empty for this particle type */
        for(i = 0; i < NumPart; i++)
        {
            if(All.MassTable[P[i].Type] != 0)
                P[i].Mass = All.MassTable[P[i].Type];
        }
    }
    
    /* zero this out, since various operations in the code will want to change particle
     masses and keeping MassTable fixed won't allow that to happen */
    for(i=0;i<6;i++) All.MassTable[i]=0;
    
    
    
#ifdef GALSF
    if(RestartFlag == 0)
    {
        if(All.MassTable[4] == 0 && All.MassTable[0] > 0)
        {
            All.MassTable[0] = 0;
            All.MassTable[4] = 0;
        }
    }
#endif
    
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE)
    if(RestartFlag == 0)
    {
        All.MassTable[2] = 0;
        All.MassTable[3] = 0;
        All.MassTable[4] = 0;
    }
#endif
    
    
    
    u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
    u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* unit conversion */
    
    if(All.InitGasTemp > 1.0e4)	/* assuming FULL ionization */
        molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
    else				/* assuming NEUTRAL GAS */
        molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);
    
    u_init /= molecular_weight;
    
    All.InitGasU = u_init;
    
    if(RestartFlag == 0)
    {
        if(All.InitGasTemp > 0)
        {
            for(i = 0; i < N_gas; i++)
            {
                if(ThisTask == 0 && i == 0)// && SphP[i].InternalEnergy == 0)
                    printf("Initializing u from InitGasTemp : InitGasTemp=%g InitGasU=%g MinEgySpec=%g SphP[0].InternalEnergy=%g\n",
                           All.InitGasTemp,All.InitGasU,All.MinEgySpec,SphP[i].InternalEnergy);
                
                SphP[i].InternalEnergy = All.InitGasU;
            }
        }
    }
    
    for(i = 0; i < N_gas; i++)
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(ThisTask == 0)
    {
        printf("reading done.\n");
        fflush(stdout);
    }
    
    if(ThisTask == 0)
    {
        printf("Total number of particles :  %d%09d\n\n",
               (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
        fflush(stdout);
    }
    
    CPU_Step[CPU_SNAPSHOT] += measure_time();
}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
    int n, k;
    MyInputFloat *fp;
    MyInputPosFloat *fp_pos;
    MyIDType *ip;
    float *fp_single;
    
    fp = (MyInputFloat *) CommBuffer;
    fp_pos = (MyInputPosFloat *) CommBuffer;
    fp_single = (float *) CommBuffer;
    ip = (MyIDType *) CommBuffer;
    
    switch (blocknr)
    {
        case IO_POS:		/* positions */
            for(n = 0; n < pc; n++)
                for(k = 0; k < 3; k++)
                {
                    P[offset + n].Pos[k] = *fp_pos++;
                    //P[offset + n].Pos[k] += 0.5*All.BoxSize; /* manually turn on for some ICs */
                }
            
            for(n = 0; n < pc; n++)
                P[offset + n].Type = type;	/* initialize type here as well */
            break;
            
        case IO_VEL:		/* velocities */
            for(n = 0; n < pc; n++)
                for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
                /* scaling v to use same IC's for different cosmologies */
                    if(RestartFlag == 0)
                        P[offset + n].Vel[k] = (*fp++) * All.VelIniScale;
                    else
                        P[offset + n].Vel[k] = *fp++;
#else
            P[offset + n].Vel[k] = *fp++;
#endif
            break;
            
        case IO_ID:		/* particle ID */
            for(n = 0; n < pc; n++)
                P[offset + n].ID = *ip++;
            break;
            
            
        case IO_CHILD_ID:		// particle child ID //
            if(RestartFlag == 2)
            {
                for(n = 0; n < pc; n++)
                    P[offset + n].ID_child_number = *ip++;
            }
            break;

        case IO_GENERATION_ID:		// particle generation ID //
            if(RestartFlag == 2)
            {
                for(n = 0; n < pc; n++)
                    P[offset + n].ID_generation = *ip++;
            }
            break;

        case IO_MASS:		/* particle mass */
            for(n = 0; n < pc; n++)
                P[offset + n].Mass = *fp++;
            break;
            
            
        case IO_SHEET_ORIENTATION:	/* initial particle sheet orientation */
            break;
            
        case IO_INIT_DENSITY:	/* initial stream density */
            
        case IO_CAUSTIC_COUNTER:	/* initial caustic counter */
            
        case IO_SECONDORDERMASS:
            for(n = 0; n < pc; n++)
            {
                P[offset + n].OldAcc = P[offset + n].Mass;	/* use this to temporarily store the masses in the 2plt IC case */
                P[offset + n].Mass = *fp++;
            }
            break;
            
        case IO_U:			/* temperature */
            for(n = 0; n < pc; n++)
                SphP[offset + n].InternalEnergy = *fp++;
            break;
            
        case IO_RHO:		/* density */
            for(n = 0; n < pc; n++)
                SphP[offset + n].Density = *fp++;
            break;
            
        case IO_NE:		/* electron abundance */
#if defined(COOLING) || defined(FLAG_NOT_IN_PUBLIC_CODE)
            for(n = 0; n < pc; n++)
                SphP[offset + n].Ne = *fp++;
#endif
            break;
            
            
        case IO_HSML:		/* gas kernel length */
            for(n = 0; n < pc; n++)
                PPP[offset + n].Hsml = *fp++;
            break;
            
        case IO_DELAYTIME:
#ifdef GALSF_FB_LUPI
	    for(n = 0; n < pc; n++)
		SphP[offset + n].DelayTimeCoolingSNe = *fp++;
#endif
            break;
            
	case IO_STELLARINITMASS:
#ifdef GALSF_FB_LUPI
            for(n = 0; n < pc; n++)
		P[offset + n].StellarInitMass = *fp++;
#endif
	    break;
	    
        case IO_AGE:		/* Age of stars */
#ifdef GALSF
            for(n = 0; n < pc; n++)
                P[offset + n].StellarAge = *fp++;
#endif
            break;
            
        case IO_GRAINSIZE:
            break;
            
        case IO_Z:			/* Gas and star metallicity */
#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(GRACKLE_OPTS)
            for(n = 0; n < pc; n++)
            {
                for(k = 0; k < NUM_METAL_SPECIES; k++)
                    P[offset + n].Metallicity[k] = *fp++;
            }

#endif
            break;
            
        case IO_VRMS:		/* Turbulence on kernel scale */
            break;
        case IO_VBULK:
            break;
        case IO_VTAN:
            break;
        case IO_VRAD:
            break;
        case IO_VDIV:
            break;
        case IO_VROT:
            break;
        case IO_VORT:
            break;
        case IO_TRUENGB:
            break;
            
        case IO_BFLD:		/* Magnetic field */
#ifdef MAGNETIC
            for(n = 0; n < pc; n++)
                for(k = 0; k < 3; k++) {SphP[offset + n].BPred[k] = *fp++;}
            SphP[offset + n].divB = 0;
#ifdef DIVBCLEANING_DEDNER
            SphP[offset + n].Phi = 0;
            SphP[offset + n].PhiPred = 0;
#endif
#endif
            break;
            
	case IO_BHSTOREDENERGY:
#ifdef BH_LUPI
	    for(n = 0; n < pc; n++)
		P[offset + n].BH_StoredEnergy = *fp++;
#endif
	    break;
            
        case IO_BHMASS:
            break;
            
        case IO_BH_DIST:
            break;
            
        case IO_BHMASSALPHA:
            break;
            
        case IO_BHMDOT:
            break;
            
        case IO_BHPROGS:
#ifdef BH_COUNTPROGS
            for(n = 0; n < pc; n++)
                P[offset + n].BH_CountProgs = *fp++;
#endif
            break;
            
        case IO_BHMBUB:
            break;
            
        case IO_BHMINI:
            break;
            
        case IO_BHMRAD:
            break;
            
        case IO_EOSTEMP:
#ifdef EOS_CARRIES_TEMPERATURE
            for(n = 0; n < pc; n++)
                SphP[offset + n].Temperature = *fp++;
#endif
            break;
            
        case IO_EOSABAR:
#ifdef EOS_CARRIES_ABAR
            for(n = 0; n < pc; n++)
                SphP[offset + n].Abar = *fp++;
#endif
            break;
            
        case IO_EOSYE:
#ifdef EOS_CARRIES_YE
            for(n = 0; n < pc; n++)
                SphP[offset + n].Ye = *fp++;
#endif
            break;
            
        case IO_RADGAMMA:
            break;
            
        case IO_DMHSML:
            break;
            
        case IO_DMDENSITY:
            break;
            
        case IO_DMVELDISP:
            break;
            
        case IO_DMHSML_V:
            break;
            
        case IO_DMDENSITY_V:
            break;
            
        case IO_CHEM:		/* Chemical abundances */
            break;
            
            
            /* adaptive softening parameters */
        case IO_AGS_SOFT:
#if defined (ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTGRAVSOFT)
            for(n = 0; n < pc; n++)
                PPP[offset + n].AGS_Hsml = *fp++;
#endif
            break;

        case IO_AGS_ZETA:
#if defined (ADAPTIVE_GRAVSOFT_FORALL) && defined(AGS_OUTPUTZETA)
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].AGS_Zeta = *fp++;
#endif
            break;

        case IO_AGS_OMEGA:
        case IO_AGS_CORR:
        case IO_AGS_NGBS:
            break;

            
            /* the other input fields (if present) are not needed to define the
             initial conditions of the code */
            
        case IO_SFR:
#ifdef GALSF
            for(n = 0; n < pc; n++)
             	PPPZ[offset + n].Sfr = *fp++;
	    break;
#endif

        case IO_HSMS:
#ifdef GALSF_FB_LUPI
            for(n = 0; n < pc; n++)
                PPP[offset + n].Hsml = *fp++;
	    break;
#endif
        case IO_POT:
        case IO_ACCEL:
        case IO_DTENTR:
        case IO_STRESSDIAG:
        case IO_STRESSOFFDIAG:
        case IO_RAD_ACCEL:
        case IO_DISTORSIONTENSORPS:
        case IO_HeHII:
        case IO_DI:
        case IO_DII:
        case IO_HD:
        case IO_HM:
        case IO_H2II:
        case IO_H2I:
        case IO_HeIII:
        case IO_HeII:
        case IO_HeI:
        case IO_HII:
        case IO_NH:
        case IO_STRESSBULK:
        case IO_SHEARCOEFF:
        case IO_TSTP:
        case IO_DBDT:
        case IO_IMF:
        case IO_COSMICRAY_ENERGY:
        case IO_DIVB:
        case IO_ABVC:
        case IO_COOLRATE:
        case IO_CONDRATE:
        case IO_DENN:
        case IO_AMDC:
        case IO_PHI:
        case IO_GRADPHI:
        case IO_TIDALTENSORPS:
        case IO_ROTB:
        case IO_FLOW_DETERMINANT:
        case IO_STREAM_DENSITY:
        case IO_PHASE_SPACE_DETERMINANT:
        case IO_ANNIHILATION_RADIATION:
        case IO_PRESSURE:
        case IO_EDDINGTON_TENSOR:
        case IO_LAST_CAUSTIC:
        case IO_ACRB:
        case IO_VSTURB_DISS:
        case IO_VSTURB_DRIVE:
        case IO_MG_PHI:
        case IO_MG_ACCEL:
	    break;
        case IO_grHI:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHI = *fp++;
	    break;
#endif
        case IO_grHII:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHII = *fp++;
	    break;
#endif
        case IO_grHM:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHM = *fp++;
	    break;
#endif
        case IO_grHeI:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHeI = *fp++;
	    break;
#endif
        case IO_grHeII:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHeII = *fp++;
	    break;
#endif
        case IO_grHeIII:
#if GRACKLE_CHEMISTRY>0
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHeIII = *fp++;
	    break;
#endif
        case IO_grH2I:
#if GRACKLE_CHEMISTRY>1
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grH2I = *fp++;
	    break;
#endif
        case IO_grH2II:
#if GRACKLE_CHEMISTRY>1
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grH2II = *fp++;
	    break;
#endif
        case IO_grDI:
#if GRACKLE_CHEMISTRY>2
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grDI = *fp++;
	    break;
#endif
        case IO_grDII:
#if GRACKLE_CHEMISTRY>2
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grDII = *fp++;
	    break;
#endif
        case IO_grHDI:
#if GRACKLE_CHEMISTRY>2
            for(n = 0; n < pc; n++)
                PPPZ[offset + n].grHDI = *fp++;
	    break;
#endif
	    break;
            
        case IO_LASTENTRY:
            endrun(220);
            break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char *fname, int readTask, int lastTask)
{
    size_t blockmaxlen;
    long long i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task, nall, nread, nstart, npart;
    int blksize1, blksize2, type, bnr, bytes_per_blockelement, nextblock, typelist[6];
    MPI_Status status;
    FILE *fd = 0;
    char label[4], buf[500];
    enum iofields blocknr;
    size_t bytes;
    
#ifdef HAVE_HDF5
    int rank, pcsum;
    hid_t hdf5_file = 0, hdf5_grp[6], hdf5_dataspace_in_file;
    hid_t hdf5_datatype = 0, hdf5_dataspace_in_memory, hdf5_dataset;
    hsize_t dims[2], count[2], start[2];
#endif
    
#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}
    
    if(ThisTask == readTask)
    {
        if(All.ICFormat == 1 || All.ICFormat == 2)
        {
            if(!(fd = fopen(fname, "r")))
            {
                printf("can't open file `%s' for reading initial conditions.\n", fname);
                endrun(123);
            }
            
            
            if(All.ICFormat == 2)
            {
                SKIP;
                my_fread(&label, sizeof(char), 4, fd);
                my_fread(&nextblock, sizeof(int), 1, fd);
                printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
                       nextblock);
                SKIP2;
            }
            
            SKIP;
            my_fread(&header, sizeof(header), 1, fd);
            SKIP2;
            
            if(blksize1 != 256 || blksize2 != 256)
            {
                printf("incorrect header format\n");
                fflush(stdout);
                endrun(890);
                /* Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total. */
            }
        }
        
        
#ifdef HAVE_HDF5
        if(All.ICFormat == 3)
        {
            read_header_attributes_in_hdf5(fname);
            hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
            for(type = 0; type < 6; type++)
            {
                if(header.npart[type] > 0)
                {
                    sprintf(buf, "/PartType%d", type);
                    hdf5_grp[type] = H5Gopen(hdf5_file, buf);
                }
            }
        }
#endif
        
        for(task = readTask + 1; task <= lastTask; task++)
        {
            MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MPI_COMM_WORLD);
        }
        
    }
    else
    {
        MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MPI_COMM_WORLD, &status);
    }
    
#ifdef INPUT_IN_DOUBLEPRECISION
    if(header.flag_doubleprecision == 0)
    {
        if(ThisTask == 0)
            printf
            ("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
        endrun(11);
    }
#else
    if(header.flag_doubleprecision)
    {
        if(ThisTask == 0)
            printf
            ("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
        endrun(10);
    }
#endif
    
    
    if(All.TotNumPart == 0)
    {
        if(header.num_files <= 1)
            for(i = 0; i < 6; i++)
            {
                header.npartTotal[i] = header.npart[i];
                header.npartTotalHighWord[i] = 0;
            }
        
        All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);
        
        for(i = 0, All.TotNumPart = 0; i < 6; i++)
        {
            All.TotNumPart += header.npartTotal[i];
            All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
        }
        
        
        for(i = 0; i < 6; i++)
            All.MassTable[i] = header.mass[i];
        
        All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));
        All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may reside on a processor */
        All.MaxPartSph = All.MaxPart; // PFH: increasing All.MaxPartSph according to this line can allow better load-balancing in some cases. however it leads to more memory problems
        // (PFH: needed to revert the change -- i.e. INCLUDE the line above: commenting it out, while it improved memory useage, causes some instability in the domain decomposition for
        //   sufficiently irregular trees. overall more stable behavior with the 'buffer', albeit at the expense of memory )
        
        
        allocate_memory();
        
        if(!(CommBuffer = mymalloc("CommBuffer", bytes = All.BufferSize * 1024 * 1024)))
        {
            printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(2);
        }
        
        if(RestartFlag >= 2)
        {
            All.Time = All.TimeBegin = header.time;
            set_cosmo_factors_for_current_time();
        }
        
    }
    
    if(ThisTask == readTask)
    {
        for(i = 0, n_in_file = 0; i < 6; i++)
            n_in_file += header.npart[i];
        
        printf("\nreading file `%s' on task=%d (contains %lld particles.)\n"
               "distributing this file to tasks %d-%d\n"
               "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
               "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
               lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
               (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
               (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
               All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
               (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
               (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
               All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
               (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
               (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
               All.MassTable[5]);
        fflush(stdout);
    }
    
    
    ntask = lastTask - readTask + 1;
    
    
    /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */
    
    for(type = 0, nall = 0; type < 6; type++)
    {
        n_in_file = header.npart[type];
        
        n_for_this_task = n_in_file / ntask;
        if((ThisTask - readTask) < (n_in_file % ntask))
            n_for_this_task++;
        
        
        if(type == 0)
        {
            if(N_gas + n_for_this_task > All.MaxPartSph)
            {
                printf("Not enough space on task=%d for SPH particles (space for %d, need at least %lld)\n",
                       ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
                fflush(stdout);
                endrun(172);
            }
        }
        
        nall += n_for_this_task;
    }
    
    if(NumPart + nall > All.MaxPart)
    {
        printf("Not enough space on task=%d (space for %d, need at least %lld)\n",
               ThisTask, All.MaxPart, NumPart + nall);
        fflush(stdout);
        endrun(173);
    }
    
    memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
    nstart = N_gas;
    
    
    
    for(bnr = 0; bnr < 1000; bnr++)
    {
        blocknr = (enum iofields) bnr;
        
        if(blocknr == IO_LASTENTRY)
            break;
        
        if(RestartFlag == 5 && blocknr > IO_MASS)	/* if we only do power spectra, we don't need to read other blocks beyond the mass */
            continue;
        
        
        if(blockpresent(blocknr))
        {
                if(RestartFlag == 0 && blocknr > IO_U && blocknr != IO_BFLD
#ifdef READ_HSML
                   && blocknr != IO_HSML
#endif
#ifdef EOS_CARRIES_TEMPERATURE
                   && blocknr != IO_EOSTEMP
#endif
#ifdef EOS_CARRIES_ABAR
                   && blocknr != IO_EOSABAR
#endif
#ifdef EOS_CARRIES_YE
                   && blocknr != IO_EOSYE
#endif
                   )
                                continue;	/* ignore all other blocks in initial conditions */
            

            if(RestartFlag == 0 && (blocknr == IO_GENERATION_ID || blocknr == IO_CHILD_ID)) continue;
#if defined(NO_CHILD_IDS_IN_ICS) || defined(ASSIGN_NEW_IDS)
            if(blocknr == IO_GENERATION_ID || blocknr == IO_CHILD_ID) continue;
#endif
            
            
            
#ifdef B_SET_IN_PARAMS
            if(RestartFlag == 0 && blocknr == IO_BFLD)
                continue;
#endif
            
            
#ifdef ADAPTIVE_GRAVSOFT_FORALL
#ifndef AGS_OUTPUTGRAVSOFT
            if(blocknr == IO_AGS_SOFT)
                continue;
#endif
#ifndef AGS_OUTPUTZETA
            if(blocknr == IO_AGS_ZETA)
                continue;
#endif
            if(blocknr == IO_AGS_OMEGA)
                continue;
            if(blocknr == IO_AGS_NGBS)
                continue;
            if(blocknr == IO_AGS_CORR)
                continue;
#endif
            
            
#ifndef GALSF //ALupi
            if(blocknr == IO_HSMS)
                continue;
#endif
            
            if(ThisTask == readTask)
            {
                get_dataset_name(blocknr, buf);
                printf("reading block %d (%s)...\n", bnr, buf);
                fflush(stdout);
            }
            
            bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);
            
            blockmaxlen = (size_t) ((All.BufferSize * 1024 * 1024) / bytes_per_blockelement);
            
            npart = get_particles_in_block(blocknr, &typelist[0]);
            
            if(npart > 0)
            {
                if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP
                   && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
                    if(ThisTask == readTask)
                    {
                        if(All.ICFormat == 2)
                        {
                            get_Tab_IO_Label(blocknr, label);
                            find_block(label, fd);
                        }
                        
                        if(All.ICFormat == 1 || All.ICFormat == 2)
                            SKIP;
                    }
                
                for(type = 0, offset = 0, nread = 0; type < 6; type++)
                {
                    n_in_file = header.npart[type];
#ifdef HAVE_HDF5
                    pcsum = 0;
#endif
                    if(typelist[type] == 0)
                    {
                        n_for_this_task = n_in_file / ntask;
                        if((ThisTask - readTask) < (n_in_file % ntask))
                            n_for_this_task++;
                        
                        offset += n_for_this_task;
                    }
                    else
                    {
                        for(task = readTask; task <= lastTask; task++)
                        {
                            n_for_this_task = n_in_file / ntask;
                            if((task - readTask) < (n_in_file % ntask))
                                n_for_this_task++;
                            
                            if(task == ThisTask)
                                if(NumPart + n_for_this_task > All.MaxPart)
                                {
                                    printf("too many particles. %d %lld %d\n", NumPart, n_for_this_task, All.MaxPart);
                                    endrun(1313);
                                }
                            
                            
                            do
                            {
                                pc = n_for_this_task;
                                
                                if(pc > (int)blockmaxlen)
                                    pc = blockmaxlen;
                                
                                if(ThisTask == readTask)
                                {
                                    if(All.ICFormat == 1 || All.ICFormat == 2)
                                    {
                                        if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY
                                           && blocknr != IO_DMVELDISP && blocknr != IO_DMHSML_V
                                           && blocknr != IO_DMDENSITY_V)
                                        {
                                            my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
                                            nread += pc;
                                        }
                                        else
                                        {
                                            nread += pc;
                                        }
                                    }
                                    
#ifdef HAVE_HDF5
                                    if(All.ICFormat == 3 && pc > 0)
                                    {
                                        get_dataset_name(blocknr, buf);
                                        hdf5_dataset = H5Dopen(hdf5_grp[type], buf);
                                        
                                        dims[0] = header.npart[type];
                                        dims[1] = get_values_per_blockelement(blocknr);
                                        if(dims[1] == 1)
                                            rank = 1;
                                        else
                                            rank = 2;
                                        
                                        hdf5_dataspace_in_file = H5Screate_simple(rank, dims, NULL);
                                        
                                        dims[0] = pc;
                                        hdf5_dataspace_in_memory = H5Screate_simple(rank, dims, NULL);
                                        
                                        start[0] = pcsum;
                                        start[1] = 0;
                                        
                                        count[0] = pc;
                                        count[1] = get_values_per_blockelement(blocknr);
                                        pcsum += pc;
                                        
                                        H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET,
                                                            start, NULL, count, NULL);
                                        
                                        switch (get_datatype_in_block(blocknr))
                                        {
                                            case 0:
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT);
                                                break;
                                            case 1:
#ifdef INPUT_IN_DOUBLEPRECISION
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                                break;
                                            case 2:
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_UINT64);
                                                break;
                                                
                                            case 3:
#if defined(INPUT_POSITIONS_IN_DOUBLE)
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                                hdf5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                                break;
                                        }
                                        
                                        H5Dread(hdf5_dataset, hdf5_datatype, hdf5_dataspace_in_memory,
                                                hdf5_dataspace_in_file, H5P_DEFAULT, CommBuffer);
                                        
                                        H5Tclose(hdf5_datatype);
                                        H5Sclose(hdf5_dataspace_in_memory);
                                        H5Sclose(hdf5_dataspace_in_file);
                                        H5Dclose(hdf5_dataset);
                                    }
#endif
                                }
                                
                                if(ThisTask == readTask && task != readTask && pc > 0)
                                    MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
                                              TAG_PDATA, MPI_COMM_WORLD);
                                
                                if(ThisTask != readTask && task == ThisTask && pc > 0)
                                    MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
                                             TAG_PDATA, MPI_COMM_WORLD, &status);
                                
                                if(ThisTask == task)
                                {
                                    empty_read_buffer(blocknr, nstart + offset, pc, type);
                                    
                                    offset += pc;
                                }
                                
                                n_for_this_task -= pc;
                            }
                            while(n_for_this_task > 0);
                        }
                    }
                }
                
                if(ThisTask == readTask)
                {
                    if(blocknr != IO_DMHSML && blocknr != IO_DMDENSITY && blocknr != IO_DMVELDISP
                       && blocknr != IO_DMHSML_V && blocknr != IO_DMDENSITY_V)
                        if(All.ICFormat == 1 || All.ICFormat == 2)
                        {
                            SKIP2;
                            
                            if(blksize1 != blksize2)
                            {
			    printf("eccomi %d\n",blocknr);
	    			fflush(stdout);

                                printf("incorrect block-sizes detected!\n");
                                printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, bnr,
                                       blksize1, blksize2);
                                if(blocknr == IO_ID)
                                {
                                    printf
                                    ("Possible mismatch of 32bit and 64bit ID's in IC file and GIZMO compilation !\n");
                                }
                                fflush(stdout);
                                endrun(1889);
                            }
                        }
                }
            }
        }
    }
    
    
    
    
    for(type = 0; type < 6; type++)
    {
        n_in_file = header.npart[type];
        
        n_for_this_task = n_in_file / ntask;
        if((ThisTask - readTask) < (n_in_file % ntask))
            n_for_this_task++;
        
        NumPart += n_for_this_task;
        
        if(type == 0)
            N_gas += n_for_this_task;
    }
    
    if(ThisTask == readTask)
    {
        if(All.ICFormat == 1 || All.ICFormat == 2)
            fclose(fd);
#ifdef HAVE_HDF5
        if(All.ICFormat == 3)
        {
            for(type = 5; type >= 0; type--)
                if(header.npart[type] > 0)
                    H5Gclose(hdf5_grp[type]);
            H5Fclose(hdf5_file);
        }
#endif
    }
    
}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
    FILE *fd;
    char buf[200], buf1[200];
    int dummy;
    
    sprintf(buf, "%s.%d", fname, 0);
    sprintf(buf1, "%s", fname);
    
    if(All.ICFormat == 3)
    {
        sprintf(buf, "%s.%d.hdf5", fname, 0);
        sprintf(buf1, "%s.hdf5", fname);
    }
    
#ifndef  HAVE_HDF5
    if(All.ICFormat == 3)
    {
        if(ThisTask == 0)
            printf("Code wasn't compiled with HDF5 support enabled!\n");
        endrun(0);
    }
#endif
    
    header.num_files = 0;
    
    if(ThisTask == 0)
    {
        if((fd = fopen(buf, "r")))
        {
            if(All.ICFormat == 1 || All.ICFormat == 2)
            {
                if(All.ICFormat == 2)
                {
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                }
                
                my_fread(&dummy, sizeof(dummy), 1, fd);
                my_fread(&header, sizeof(header), 1, fd);
                my_fread(&dummy, sizeof(dummy), 1, fd);
            }
            fclose(fd);
            
#ifdef HAVE_HDF5
            if(All.ICFormat == 3)
                read_header_attributes_in_hdf5(buf);
#endif
        }
    }
    
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    if(header.num_files > 0)
        return header.num_files;
    
    if(ThisTask == 0)
    {
        if((fd = fopen(buf1, "r")))
        {
            if(All.ICFormat == 1 || All.ICFormat == 2)
            {
                if(All.ICFormat == 2)
                {
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                    my_fread(&dummy, sizeof(dummy), 1, fd);
                }
                
                my_fread(&dummy, sizeof(dummy), 1, fd);
                my_fread(&header, sizeof(header), 1, fd);
                my_fread(&dummy, sizeof(dummy), 1, fd);
            }
            fclose(fd);
            
#ifdef HAVE_HDF5
            if(All.ICFormat == 3)
                read_header_attributes_in_hdf5(buf1);
#endif
            
            header.num_files = 1;
        }
    }
    
    MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    if(header.num_files > 0)
        return header.num_files;
    
    if(ThisTask == 0)
    {
        printf("\nCan't find initial conditions file.");
        printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
        fflush(stdout);
    }
    
    endrun(0);
    return 0;
}

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(FLAG_NOT_IN_PUBLIC_CODE_RESHUFFLE_CATALOGUE)
void get_particle_numbers(char *fname, int num_files)
{
    char buf[1000];
    int blksize1, blksize2;
    char label[4];
    int nextblock;
    int i, j;
    
    printf("num_files=%d\n", num_files);
    
    for(i = 0; i < num_files; i++)
    {
        if(num_files > 1)
        {
            sprintf(buf, "%s.%d", fname, i);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.%d.hdf5", fname, i);
        }
        else
        {
            sprintf(buf, "%s", fname);
            if(All.ICFormat == 3)
                sprintf(buf, "%s.hdf5", fname);
        }
        
#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}
        
        if(All.ICFormat == 1 || All.ICFormat == 2)
        {
            FILE *fd;
            
            if(!(fd = fopen(buf, "r")))
            {
                printf("can't open file `%s' for reading initial conditions.\n", buf);
                endrun(1239);
            }
            
            if(All.ICFormat == 2)
            {
                SKIP;
                my_fread(&label, sizeof(char), 4, fd);
                my_fread(&nextblock, sizeof(int), 1, fd);
                SKIP2;
            }
            
            SKIP;
            my_fread(&header, sizeof(header), 1, fd);
            SKIP2;
            if(blksize1 != 256 || blksize2 != 256)
            {
                printf("incorrect header format\n");
                fflush(stdout);
                endrun(890);
            }
            fclose(fd);
        }
        
#ifdef HAVE_HDF5
        if(All.ICFormat == 3)
        {
            read_header_attributes_in_hdf5(buf);
        }
#endif
        
        NumPartPerFile[i] = 0;
        
        for(j = 0; j < 6; j++)
        {
                NumPartPerFile[i] += header.npart[j];
        }
        
        printf("File=%4d:  NumPart= %d\n", i, (int) (NumPartPerFile[i]));
    }
    
    
    long long n, sum;
    
    for(i = 0, sum = 0; i < num_files; i++)
    {
        n = NumPartPerFile[i];
        
        NumPartPerFile[i] = sum;
        
        sum += n;
    }
}
#endif




/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
                     int *last)
{
    int ntask, filesleft, filesright, tasksleft;
    
    if(nfiles > 1)
    {
        ntask = lasttask - firsttask + 1;
        
        filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
        if(filesleft <= 0)
            filesleft = 1;
        if(filesleft >= nfiles)
            filesleft = nfiles - 1;
        
        filesright = nfiles - filesleft;
        
        tasksleft = ntask / 2;
        
        distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
        distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
                        last);
    }
    else
    {
        if(ThisTask >= firsttask && ThisTask <= lasttask)
        {
            *filenr = firstfile;
            *master = firsttask;
            *last = lasttask;
        }
    }
}



#ifdef HAVE_HDF5
void read_header_attributes_in_hdf5(char *fname)
{
    hid_t hdf5_file, hdf5_headergrp, hdf5_attribute;
    
    hdf5_file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    hdf5_headergrp = H5Gopen(hdf5_file, "/Header");
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_ThisFile");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, header.npart);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPart_Total_HighWord");
    H5Aread(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "MassTable");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
    H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumFilesPerSnapshot");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.num_files);
    H5Aclose(hdf5_attribute);
    
    hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Flag_DoublePrecision");
    H5Aread(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision);
    H5Aclose(hdf5_attribute);
    
    H5Gclose(hdf5_headergrp);
    H5Fclose(hdf5_file);
}
#endif






/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void find_block(char *label, FILE * fd)
{
    unsigned int blocksize = 0, blksize;
    char blocklabel[5] = { "    " };
    
#define FBSKIP  {my_fread(&blksize,sizeof(int),1,fd);}
    
    rewind(fd);
    
    while(!feof(fd) && blocksize == 0)
    {
        FBSKIP;
        if(blksize != 8)
        {
            printf("Incorrect Format (blksize=%u)!\n", blksize);
            exit(1891);
        }
        else
        {
            my_fread(blocklabel, 4 * sizeof(char), 1, fd);
            my_fread(&blocksize, sizeof(int), 1, fd);
            /*
             printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
             label[0],label[1],label[2],label[3],blocklabel,blocksize);
             */
            FBSKIP;
            if(strncmp(label, blocklabel, 4) != 0)
            {
                fseek(fd, blocksize, 1);
                blocksize = 0;
            }
        }
    }
    if(feof(fd))
    {
        printf("Block '%c%c%c%c' not found !\n", label[0], label[1], label[2], label[3]);
        fflush(stdout);
        endrun(1890);
    }
}
