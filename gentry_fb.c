#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"
#include "kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif


#ifdef WINDS
#define NWT (2)
#else
#define NWT (1)
#endif


/* Routines for thermal feedback by SNe as in Gentry et al. 2016 */

/*
 * This file was written by Eric Gentry (egentry@ucsc.edu) for GIZMO.
 * This was based on Alessandro's creasey.c file
 */

// to do: remove the assumption that SNe come from star particles which
// can move and then inject energy within their own particles

void neighbour_loop(int,int, double*);


void gentry_fb_calc()
{

  // double *dtWinds;
  int i,j;
  int first_SN = All.N_SNe;
  int last_SN  = -1;
  int first_OB = All.N_SNe;
  int last_OB  = -1;

  if(All.TimeStep == 0) return;
  if(ThisTask==0) printf("Start SN/OB release\n");
  fflush(stdout);

  // dtWinds =(double*)mymalloc("dtWinds",All.N_SNe*sizeof(double));
    
  // // Find active stars and SNe
  for(i=0;i<All.N_SNe;i++)
  {

#ifdef WINDS
    //OB stars
    if(All.Time<All.SN_time[i])
    {
      if(i<=first_OB) first_OB=i;
      if(i>=last_OB)   last_OB=i;
      // dtWinds[i]=All.TimeStep;
      // if(All.Time>=All.SN_time[i] && All.Time-All.TimeStep<All.SN_time[i])
      //   dtWinds[i]=DMIN(All.TimeStep,All.SN_time[i]-All.Time+All.TimeStep);
    }
#endif


    //SNe
    if(All.Time>All.SN_time[i] && All.Time-All.TimeStep<All.SN_time[i])
    {
      if(i<=first_SN) first_SN=i;
      if(i>=last_SN)   last_SN=i;
    }
  }

  if(ThisTask==0)
    printf("Stars: %d %d %d %d\n",first_OB,last_OB,first_SN,last_SN);
  fflush(stdout);

  double my_weight[NWT],total_weight[NWT];
  for(j=0;j<NWT;j++)
    my_weight[j]=total_weight[j]=0.;
  for(i=0; i<All.N_SNe; i++)
  {
    // calculate myweight
    if(i>=IMIN(first_OB,first_SN) && i<=IMAX(last_OB,last_SN)) neighbour_loop(0,i,my_weight);
    //printf("%d %g\n",ThisTask,myweight);
    //fflush(stdout);

    MPI_Allreduce(&my_weight,&total_weight,NWT,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //printf("MY error: %g\n",totWeight);
    if(i>=first_SN && i<=last_SN) neighbour_loop(2,i,total_weight);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // myfree(dtWinds);
}


/* neighbour_loop -- Calculates weight or applies feedback for neighbourings

    Mode 0:
      Calculates my_weight for the ith massive star, for.
        - if including winds, my_weight is an array of 2 values
          - my_weight[0] : weight for the SN feedback
          - my_weight[1] : weight for the stellar wind feedback
        - if not including winds, it is the same as above,
          but only an array of 1 value, my_weight[0] (SN feedback)

        - Note, this only includes neighbours on the local processor
          To get the total weight across the domain, need an MPI_ALLreduce

      Returns the weight through double *out_weight

    Mode 1:
      Add winds

    Mode 2: 
      Add SNe

*/
void neighbour_loop(int mode,int iSN,double *out_weight)
{
    int j, k, n;

    double rSN = 2.0 * 3.085678e18 / All.UnitLength_in_cm; // injection radius
#ifdef WINDS
    double pdotWIND = 920*1.989e38 / All.UnitMass_in_g 
                    / All.UnitVelocity_in_cm_per_s * All.UnitTime_in_Megayears;
    double mdotWIND = 0.1*1.989e33/ All.UnitMass_in_g*All.UnitTime_in_Megayears;
    double edotWIND = 1.0e50 / All.UnitEnergy_in_cgs*All.UnitTime_in_Megayears;
#endif

    double eSN = 1.0e51 / All.UnitEnergy_in_cgs;
    double mSN = All.SN_mass[iSN];
#ifdef GRACKLE_OPTS
    double zSN = All.SN_mass_Z[iSN];
#endif

    double e_shock   = 0;   // amount of internal energy to inject
    double p_shock   = 0;   // amount of momentum to inject
    double m_shock   = 0;   // amount of mass to inject
#ifdef GRACKLE_OPTS
    double met_shock = 0;   // amount of metals to inject
#endif

    /* Load the data for the particle injecting feedback */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
  
    double pos[3];
    pos[0] = All.SN_position_x[iSN];
    pos[1] = All.SN_position_y[iSN];
    pos[2] = All.SN_position_z[iSN];

    double my_weight[NWT];
    for(k=0; k<NWT; k++)
      my_weight[k] = 0.;
    if(mode>0)
      for(k=0; k<NWT; k++)
        my_weight[k] = out_weight[k];

    if(ThisTask == 0)
    {
      // if(mode==1)
      // printf("OB star releasing wind at %g %g %g\n",pos[0],pos[1],pos[2]);
      // else
      if(mode==2)
        printf("SN exploding at %g %g %g at time %g (supposed to explode at %g)\n",
               pos[0],pos[1],pos[2],All.Time,All.SN_time[iSN]);
    }
    fflush(stdout);

    for(j = 0; j<NTask; j++)
    {
      Send_count[j] = 0;
      Exportflag[j] = -1;
    }

    // double total_mom=0;


    int startnode = All.MaxPart;
    while(startnode >= 0)
    {
      int dummy;
      int numngb_inbox = ngb_treefind_variable(pos,10*rSN,-1,&startnode,
                                               0,&dummy,&dummy);
          
      if(numngb_inbox < 0)
      {
        printf("I reached %g and I stopped\n",my_weight[0]);
        fflush(stdout);
        return;
      }
      for(n = 0; n < numngb_inbox; n++)
      {
        j = Ngblist[n];
        if(P[j].Type != 0) continue; // require a gas particle //
        if(P[j].Mass <= 0) continue; // require the particle has mass //
         
        double d_pos[3]; // vector distance from star center      
        for(k=0; k<3; k++) {d_pos[k] = pos[k] - P[j].Pos[k];}
#ifdef PERIODIC     /* find the closest image in the given box size  */
        NEAREST_XYZ(d_pos[0],d_pos[1],d_pos[2],1);
#endif
        double r2=0; 
        for(k=0; k<3; k++) {r2 += d_pos[k]*d_pos[k];}
        if(r2<=0) continue; // same particle //
        if(r2 >= 100*rSN*rSN) continue; // outside kernel //
  
        // compute kernel weight: gaussian in radius, 
        //   proportional to volume / number of neighbours
        double wk = 1.0/pow(2*M_PI*pow(rSN,2),1.5)*exp(-r2/(2*rSN*rSN))
                    *pow(Get_Particle_Size(j),3);
#ifdef WINDS
        // compute kernel weight for momentum. Uglier functional form.
        double wk_p[3];
        for(k=0; k<3; k++)
          wk_p[k] = 1.0/(8*M_PI*pow(rSN,4))*exp(-r2/(2*rSN*rSN))
                    *pow(Get_Particle_Size(j),3)*d_pos[k];
#endif
        if(wk==0) continue; //useless particle
  
        //Evaluate neighbour number
        if(mode == 0)
        {
          my_weight[0] += wk;
#ifdef WINDS
          my_weight[1] += sqrt(wk_p[0]*wk_p[0] 
                             + wk_p[1]*wk_p[1] 
                             + wk_p[2]*wk_p[2]);
#endif
          continue;
        }
        // if(total_weight[i]==0) printf("Error: total_weight[0]==0 \n");
        if(my_weight[0]==0) printf("Error: my_weight[0]==0 \n");
        fflush(stdout);
        // endrun(4321);
        wk /= my_weight[0];
#ifdef WINDS
        if(mode==1) //Winds
        {
          e_shock   = 0;//edotWIND*All.TimeStep;
          p_shock   = pdotWIND * All.TimeStep;
          m_shock   = mdotWIND * All.TimeStep;
#ifdef GRACKLE_OPTS
          met_shock = P[iSN].Metallicity[0] * m_shock;
#endif
        }
        else //SNe
#endif
        {
          e_shock   = eSN*wk;
          p_shock   = 0;
          m_shock   = mSN*wk;
#ifdef GRACKLE_OPTS
          met_shock = zSN*wk;
#endif
        }


        // printf("%g %g -> %g\n",All.SN_position_x[iSN],P[j].Pos[0],wk_p[0]);
        // double x=P[j].Mass/(P[j].Mass+m_shock*wk);
        // double v=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]
                     // +P[j].Vel[2]*P[j].Vel[2]);
// #ifdef WINDS
//         for(k=0;k<3;k++)
//         {
//           wk_p[k]/=my_weight[1];
//           // wk_p[k+3]/=my_weight[k+4];
//           // Momentum
//           P[j].Vel[k]        =x*(P[j].Vel[k]        - p_shock*(wk_p[k])/P[j].Mass);
//           SphP[j].VelPred[k] =x*(SphP[j].VelPred[k] - p_shock*(wk_p[k])/P[j].Mass);
//         }

//         // double vnew=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]
//                         // +P[j].Vel[2]*P[j].Vel[2]);
//         // printf("Delta v: %g / Delta p: %g of %g\n",
//                // vnew-v,P[j].Mass*(vnew-v),p_shock);
//         // fflush(stdout);
//         // total_mom+=p_shock*sqrt(pow(wk_p[0],2)+pow(wk_p[1],2)+pow(wk_p[2],2));
// #endif

        // modelling this section after `merge_particles_ij()`
        // If you naively change P[j].Mass and SphP[j].InternalEnergy,
        // then the amount of energy you add will actually depend on the
        // timestep which is determined *before* the energy is added.
        // 
        // That is bad.

        // i is the ejecta; j is the existing cell

        double mtot = P[j].Mass + m_shock;
        double wt_i = m_shock   / mtot;
        double wt_j = P[j].Mass / mtot;

        double dm_j=0,de_j=0,dp_j[3],dm_ij,de_ij,dp_ij[3];
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        dm_j += SphP[j].DtMass;
#endif
        de_j = P[j].Mass * SphP[j].DtInternalEnergy + dm_j*SphP[j].InternalEnergy;
        for(k=0;k<3;k++)
        {
            dp_j[k] = P[j].Mass * SphP[j].HydroAccel[k] \
                      + dm_j * SphP[j].VelPred[k] / All.cf_atime;
            de_j += dp_j[k] * SphP[j].VelPred[k] / All.cf_atime 
                    - 0.5 * dm_j * SphP[j].VelPred[k] * SphP[j].VelPred[k] * All.cf_a2inv;
            dp_ij[k] = dp_j[k];
        }
        dm_ij = dm_j;
        de_ij = de_j;

#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        SphP[j].MassTrue += m_shock;
        SphP[j].DtMass=dm_ij;
        SphP[j].dMass = SphP[j].dMass; // do I have to add in new mass here?
#endif

        /* make sure to update the conserved variables correctly: 
           mass and momentum are easy, energy is non-trivial */
        // double egy_old = 0;
        // internal energy //
        // egy_old += mtot * (wt_j*SphP[j].InternalEnergy + wt_i*e_shock/m_shock);
        double pos_new_xyz[3], dp[3];
        /* for periodic boxes, we need to (arbitrarily) pick one position as our coordinate center. we pick i. then everything defined in 
            position differences relative to i. the final position will be appropriately box-wrapped after these operations are completed */
        for(k=0;k<3;k++) {dp[k]=P[j].Pos[k]-pos[k];}
#ifdef PERIODIC
        NEAREST_XYZ(dp[0],dp[1],dp[2],-1);
#endif
        for(k=0;k<3;k++) {pos_new_xyz[k] = pos[k] + wt_j * dp[k];}

        // for(k=0;k<3;k++)
        // {
            // egy_old += mtot*wt_j * 0.5 * P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv; // kinetic energy (j) //
            // egy_old += mtot*wt_i * 0.5 * 0; // kinetic energy (i) //
            // // gravitational energy terms need to be added (including work for moving particles 'together') //
            // // Egrav = m*g*h = m * (-grav_acc) * (position relative to zero point) //
            // egy_old += mtot*wt_j * (pos[k]+dp[k] - pos_new_xyz[k])*All.cf_atime 
            //             * (-P[j].GravAccel[k])*All.cf_a2inv; // work (j) //
            // egy_old += mtot*wt_i * (pos[k]       - pos_new_xyz[k])*All.cf_atime 
            //             * (-0)*All.cf_a2inv;                 // work (i) //
// #ifdef HYDRO_MESHLESS_FINITE_VOLUME
//             SphP[j].GravWorkTerm[k] = 0; // since we're accounting for the work above and dont want to accidentally double-count //
// #endif
        // }
        SphP[j].InternalEnergy     =   wt_j*SphP[j].InternalEnergy 
                                     + wt_i*(e_shock / m_shock);
        
        // I'm not so sure about this step
        SphP[j].InternalEnergyPred =  wt_j*SphP[j].InternalEnergyPred 
                                    + wt_i*(e_shock / m_shock);
        double p_old_i[3],p_old_j[3];
        for(k=0;k<3;k++)
        {
            p_old_i[k] = 0;
            p_old_j[k] = P[j].Mass * P[j].Vel[k];
        }
        for(k=0;k<3;k++)
        {
            P[j].Pos[k] = pos_new_xyz[k]; // center-of-mass conserving //
            P[j].Vel[k] = wt_j*P[j].Vel[k] + 0; // momentum-conserving //
            SphP[j].VelPred[k] = wt_j*SphP[j].VelPred[k] + 0; // momentum-conserving //
            P[j].GravAccel[k] = wt_j*P[j].GravAccel[k] + 0; // force-conserving //
        }

        /* correct our 'guess' for the internal energy with the residual 
           from exact energy conservation */
        // double egy_new = mtot * SphP[j].InternalEnergy;
        // for(k=0;k<3;k++) {egy_new += mtot * 0.5*P[j].Vel[k]*P[j].Vel[k]*All.cf_a2inv;}
        // egy_new = (egy_old - egy_new) / mtot; /* this residual needs to be put into the thermal energy */
        // if(egy_new < -0.5*SphP[j].InternalEnergy) egy_new = -0.5 * SphP[j].InternalEnergy;
        //SphP[j].InternalEnergy += egy_new; SphP[j].InternalEnergyPred += egy_new;//test during splits
        if(SphP[j].InternalEnergyPred<0.5*SphP[j].InternalEnergy) SphP[j].InternalEnergyPred=0.5*SphP[j].InternalEnergy;
        
        // now use the conserved variables to correct the derivatives to primitive variables //
        de_ij -= dm_ij * SphP[j].InternalEnergyPred;
        for(k=0;k<3;k++)
        {
            SphP[j].HydroAccel[k] = (dp_ij[k] - dm_ij * SphP[j].VelPred[k]/All.cf_atime) / mtot;
            de_ij -= mtot * SphP[j].VelPred[k]/All.cf_atime * SphP[j].HydroAccel[k] + 0.5 * dm_ij * SphP[j].VelPred[k]*SphP[j].VelPred[k]*All.cf_a2inv;
        }
        SphP[j].DtInternalEnergy = de_ij;
        // to be conservative adopt the maximum signal velocity and kernel length //
        double ejecta_sound_speed = sqrt((GAMMA/(GAMMA-1)) * e_shock/m_shock);
        SphP[j].MaxSignalVel = sqrt(SphP[j].MaxSignalVel*SphP[j].MaxSignalVel 
                                    + ejecta_sound_speed); /* need to be conservative */
        // PPP[j].Hsml = pow(pow(PPP[j].Hsml,NUMDIMS)+pow(PPP[i].Hsml,NUMDIMS),1.0/NUMDIMS); /* sum the volume of the two particles */
        // SphP[j].ConditionNumber = SphP[j].ConditionNumber + SphP[i].ConditionNumber; /* sum to be conservative */
// #ifdef ENERGY_ENTROPY_SWITCH_IS_ACTIVE
        // SphP[j].MaxKineticEnergyNgb = DMAX(SphP[j].MaxKineticEnergyNgb,SphP[i].MaxKineticEnergyNgb); /* for the entropy/energy switch condition */
// #endif

        // below, we need to take care of additional physics //
#if defined(GRACKLE_OPTS)
        /* metal-mass conserving */
        P[j].Metallicity[0] = wt_j*P[j].Metallicity[0] + wt_i*met_shock/m_shock;
#endif

        P[j].Mass = mtot;
        for(k=0;k<3;k++)
        {
            /* momentum shift for passing to tree (so we know how to move it) */
            P[j].dp[k] += P[j].Mass*P[j].Vel[k] - p_old_j[k];
        }

        SphP[j].Pressure = get_pressure(j);


        // below this is code from before


        // SphP[j].InternalEnergy = (   (SphP[j].InternalEnergy*P[j].Mass)
        //                            + (e_shock*wk) )
        //                          / (P[j].Mass + m_shock*wk);

// #ifdef GRACKLE_OPTS
//         P[j].Metallicity[0]=x*(P[j].Metallicity[0] + met_shock*wk/P[j].Mass);
// #endif
        if(P[j].Mass<m_shock*wk)
        {  
          printf("Mass: %d %g %d %g %g %g %g\n",
                 iSN,m_shock,numngb_inbox,my_weight[0],P[j].Mass,
                 sqrt(r2),wk);
          fflush(stdout);
          // // endrun(1234);
        }
        // P[j].Mass+=m_shock*wk;

      } // for(n = 0; n < numngb; n++)
      // if(iSN==7 && my_weight>0) printf("particle %d: %g\n",ThisTask,my_weight);   
    } // while(startnode >= 0)

    // if(iSN==7 && my_weight>0)  printf("particle out%d: %g %d\n",
                                      // iSN,my_weight,ThisTask);    
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("My tot: %d %g %d\n",mode,my_weight,ThisTask);
    // printf("My (%d) tot momentum for %d is %g of %g\n",
            // ThisTask,iSN, total_mom,p_shock);
    //fflush(stdout);

    for(k=0;k<NWT;k++)
      out_weight[k]=my_weight[k];

    myfree(Ngblist);
} // void neighbour_loop



int count_lines_in_file(char * filename)
{
    int number_of_lines = 0;

    FILE * file = fopen(filename, "r");

    char current_char;
    char previous_char = '\n';
    while ((current_char = fgetc(file)) != EOF)
    {
      if (current_char == '\n')
      {
        if (current_char != previous_char)
        {
          ++number_of_lines;
        }
      }
      previous_char = current_char;
    }
    if (previous_char != '\n')
    {
      ++number_of_lines;
    }

    fclose(file);

    return number_of_lines;
}


void read_SNe_file(char * filename, int N_SNe, 
  double *SN_time, double *SN_mass, double *SN_mass_Z, double *wind_mass)
{
    if(N_SNe <= 0) return;

    FILE * file = fopen(filename, "r");

    double _SN_time;
    double _initial_mass;
    double _SN_mass;
    double _SN_mass_Z;
    double _wind_mass;

    char tmp[1024];
    fgets(tmp, sizeof(tmp), file); // header line

    int i;
    for(i=0; i<N_SNe; ++i)
    {
        fscanf(file,"%le %le %le %le %le\n",
                &_SN_time, &_initial_mass, 
                &_SN_mass, &_SN_mass_Z,
                &_wind_mass );

        // the original file was in *reverse* time order
        // (because it was sorted by stellar initial mass)
        SN_time[  (N_SNe-1) - i]   = _SN_time   / All.UnitTime_in_s;
        SN_mass[  (N_SNe-1) - i]   = _SN_mass   / All.UnitMass_in_g;
        SN_mass_Z[(N_SNe-1) - i]   = _SN_mass_Z / All.UnitMass_in_g;
        wind_mass[(N_SNe-1) - i]   = _wind_mass / All.UnitMass_in_g;
    }

    fclose(file);

    double first_SN_time = SN_time[0];
    for(i=0; i<N_SNe; ++i)
    {
      SN_time[i] -= first_SN_time;
    }
    SN_time[0] = 3e10 / All.UnitTime_in_s; // don't add the first SN at exactly t=0;
}
