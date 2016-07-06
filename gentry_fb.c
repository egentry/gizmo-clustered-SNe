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
  int first_SN=All.N_SNe;
  int last_SN=-1;
  int first_OB=All.N_SNe;
  int last_OB=-1;

  if(All.TimeStep == 0) return;
  if(ThisTask==0) printf("Start SN/OB release\n");
  fflush(stdout);

  // dtWinds =(double*)mymalloc("dtWinds",All.N_SNe*sizeof(double));
    
  // // Find active stars and SNe
  for(i=0;i<All.N_SNe;i++)
  {
    //OB stars
    if(All.Time>All.SN_time[i])
    {
      if(i<=first_OB) first_OB=i;
      if(i>=last_OB)   last_OB=i;
      // dtWinds[i]=All.TimeStep;
      // if(All.Time>=All.SN_time[i] && All.Time-All.TimeStep<All.SN_time[i])
      //   dtWinds[i]=DMIN(All.TimeStep,All.SN_time[i]-All.Time+All.TimeStep);
    }

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
    if(i>=first_OB && i<=DMAX(last_OB,last_SN)) neighbour_loop(0,i,my_weight);
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
    double zSN = All.SN_mass_Z[iSN];

    double e_shock   = 0;   // amount of internal energy to inject
    double p_shock   = 0;   // amount of momentum to inject
    double m_shock   = 0;   // amount of mass to inject
    double met_shock = 0;   // amount of metals to inject

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
        printf("SN exploding at %g %g %g at time %g\n",
               pos[0],pos[1],pos[2],All.SN_time[iSN]);
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
         
        double dp[3]; // vector distance from star center      
        for(k=0; k<3; k++) {dp[k] = pos[k] - P[j].Pos[k];}
#ifdef PERIODIC     /* find the closest image in the given box size  */
        NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
        double r2=0; 
        for(k=0; k<3; k++) {r2 += dp[k]*dp[k];}
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
                    *pow(Get_Particle_Size(j),3)*dp[k];
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
          met_shock = P[iSN].Metallicity[0] * m_shock;
        }
        else //SNe
#endif
        {
          e_shock   = eSN;
          p_shock   = 0;
          m_shock   = mSN;
          met_shock = zSN;
        }

        // printf("%g %g %g\n",mdotWIND,m_shock,met_shock);
        SphP[j].InternalEnergy += e_shock*wk/P[j].Mass;

#ifdef WINDS
        // printf("%g %g -> %g\n",All.SN_position_x[iSN],P[j].Pos[0],wk_p[0]);
        double x=P[j].Mass/(P[j].Mass+m_shock*wk);
        // double v=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]
                     // +P[j].Vel[2]*P[j].Vel[2]);
        for(k=0;k<3;k++)
        {
          wk_p[k]/=my_weight[1];
          // wk_p[k+3]/=my_weight[k+4];
          // Momentum
          P[j].Vel[k]        =x*(P[j].Vel[k]        - p_shock*(wk_p[k])/P[j].Mass);
          SphP[j].VelPred[k] =x*(SphP[j].VelPred[k] - p_shock*(wk_p[k])/P[j].Mass);
        }
        // double vnew=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]
                        // +P[j].Vel[2]*P[j].Vel[2]);
        // printf("Delta v: %g / Delta p: %g of %g\n",
               // vnew-v,P[j].Mass*(vnew-v),p_shock);
        // fflush(stdout);
        // total_mom+=p_shock*sqrt(pow(wk_p[0],2)+pow(wk_p[1],2)+pow(wk_p[2],2));
    
        P[j].Metallicity[0]=x*(P[j].Metallicity[0]+met_shock*wk/P[j].Mass);
        if(P[j].Mass<m_shock*wk)
        {  
          printf("Mass: %d %g %d %g %g %g %g %g\n",
                 iSN,m_shock,numngb_inbox,my_weight[0],P[j].Mass,
                 sqrt(r2),wk,wk_p[0]);
          fflush(stdout);
          // // endrun(1234);
        }
        P[j].Mass+=m_shock*wk;
#endif
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
