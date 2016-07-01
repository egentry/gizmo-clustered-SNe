#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif


#ifdef WINDS
#define NWT (2)
#else
#define NWT (1)
#endif


/* Routines for thermal feedback by SNae as in Creasey et al. 2012*/

/*
 * This file was written by Alessandro Lupi (alessandro.lupi@uninsubria.it) for GIZMO.
 */
void neighbour_loop(int,int, double*);
//void neighbour_loop(int,int, double*,double*);


void creasey_calc()
{
    double *dtWinds;
    int i,j;
    //All.NSNe=8;
    int fSN=All.NSNe,lSN=-1;
    int fOB=All.NSNe,lOB=-1;

    if(All.TimeStep == 0) return;
    if(ThisTask==0) printf("Start SN/OB release\n");
    fflush(stdout);

    dtWinds =(double*)mymalloc("dtWinds",All.NSNe*sizeof(double));
    
    //First check: active stars and SNe
    for(i=0;i<All.NSNe;i++)
    {
	//OB stars
	if(All.Time>All.tSNe[i]-All.delay && All.Time-All.TimeStep<All.tSNe[i])
	{
	  if(i<=fOB) fOB=i;
	  if(i>=lOB) lOB=i;
	  dtWinds[i]=All.TimeStep;
	  if(All.Time-All.tSNe[i]+All.delay>=0)
	    dtWinds[i]=DMIN(All.TimeStep,All.Time-All.tSNe[i]+All.delay);
          if(All.Time>=All.tSNe[i] && All.Time-All.TimeStep<All.tSNe[i])
            dtWinds[i]=DMIN(All.TimeStep,All.tSNe[i]-All.Time+All.TimeStep);
	}

	//SNe
	if(All.Time>All.tSNe[i] && All.Time-All.TimeStep<All.tSNe[i])
	{
	  if(i<=fSN) fSN=i;
	  if(i>=lSN) lSN=i;
	}
    }
    if(ThisTask==0)
	printf("Stars: %d %d %d %d\n",fOB,lOB,fSN,lSN);
    fflush(stdout);
    //getchar();

    double myweight[NWT],totWeight[NWT];
    for(j=0;j<NWT;j++)
	myweight[j]=totWeight[j]=0.;
    for(i=0;i<All.NSNe;i++)
    {
	if(i>=fOB && i<=DMAX(lOB,lSN)) neighbour_loop(0,i,myweight);
	//printf("%d %g\n",ThisTask,myweight);
	//fflush(stdout);
	MPI_Allreduce(&myweight,&totWeight,NWT,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	//printf("MY error: %g\n",totWeight);
	//getchar();
	//if(i>=fOB && i<=lOB) neighbour_loop(1,i,totWeight);
	if(i>=fSN && i<=lSN) neighbour_loop(2,i,totWeight);
	MPI_Barrier(MPI_COMM_WORLD);
    }

    myfree(dtWinds);
}





void neighbour_loop(int mode,int iSN,double *outw)
{
    int startnode, numngb_inbox, listindex = 0,dummy;
    int j, k, n;
    double r2,mywt[NWT];
    double wk;
#ifdef WINDS
    double wk_p[3];
#endif
    double dp[3];

    double pos[3];
    double rSN = 2.0 * 3.085678e18 / All.UnitLength_in_cm;
#ifdef WINDS
    double pdotWIND = 920*1.989e38 / All.UnitMass_in_g / All.UnitVelocity_in_cm_per_s * All.UnitTime_in_Megayears;
    double mdotWIND = 0.1*1.989e33/ All.UnitMass_in_g*All.UnitTime_in_Megayears;
    double edotWIND = 1.0e50 / All.UnitEnergy_in_cgs*All.UnitTime_in_Megayears;
#endif
    double zfSN = 0.2;
    double eSN = 1.0e51 / All.UnitEnergy_in_cgs;
    double mSN = 10.0*1.989e33/ All.UnitMass_in_g;
    double zSN = zfSN*mSN*1.989e33/ All.UnitMass_in_g;

    double e_shock=0;
    double p_shock=0;
    double m_shock=0;
    double met_shock=0;
    /* Load the data for the particle injecting feedback */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
//	printf("%d %d\n",fstar,lstar);
    //fstar=lstar=8;
	
    pos[0]=All.xSNe[iSN];
    pos[1]=All.ySNe[iSN];
    pos[2]=All.zSNe[iSN];
    for(k=0;k<NWT;k++)
    	mywt[k]=0.;
    if(mode>0)
        for(k=0;k<NWT;k++)
		mywt[k]=outw[k];

    if(ThisTask == 0)
	{
	 //if(mode==1)
	 //printf("OB star releasing wind at %g %g %g\n",pos[0],pos[1],pos[2]);
	 //else
	 if(mode==2)printf("SNa exploding at %g %g %g at time %g\n",pos[0],pos[1],pos[2],All.tSNe[iSN]);
	}
    fflush(stdout);
    //double totw=0,totwe=0;
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }

    double tot_mom=0;

    startnode = All.MaxPart;
    while(startnode >= 0)
        {
           while(startnode >= 0)
           {
              numngb_inbox = ngb_treefind_variable(pos,10*rSN,-1,&startnode,0,&dummy,&dummy);
            
              if(numngb_inbox < 0)
		{
		  printf("I reached %g and I stopped\n",mywt[0]);
		  fflush(stdout);
                  return;
		}
              for(n = 0; n < numngb_inbox; n++)
              {
                j = Ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                
                for(k=0; k<3; k++) {dp[k] = pos[k] - P[j].Pos[k];}
#ifdef PERIODIC			/* find the closest image in the given box size  */
                NEAREST_XYZ(dp[0],dp[1],dp[2],1);
#endif
                r2=0; for(k=0;k<3;k++) {r2 += dp[k]*dp[k];}
                if(r2<=0) continue; // same particle //
                

		if(r2 >= 100*rSN*rSN) continue; // outside kernel //
		
		wk = 1.0/pow(2*M_PI*pow(rSN,2),1.5)*exp(-r2/(2*rSN*rSN))*pow(Get_Particle_Size(j),3);
		//if(wk<1.0e-10) continue;
#ifdef WINDS
		for(k=0;k<3;k++)
		 wk_p[k] = 1.0/(8*M_PI*pow(rSN,4))*exp(-r2/(2*rSN*rSN))*pow(Get_Particle_Size(j),3)*dp[k];

#endif
		if(wk==0) continue; //useless particle
		
	        //Evaluate neighbour number
		if(mode == 0)
		{
			mywt[0] += wk;
			//printf("Valutazione: %g %g\n",wk,totWt[iSN]);
#ifdef WINDS
			//for(k=0;k<3;k++)
			// mywt[k+1]+= wk_p[k];
			mywt[1] += sqrt(wk_p[0]*wk_p[0]+wk_p[1]*wk_p[1]+wk_p[2]*wk_p[2]);
#endif
		        //printf("Wk: %g\n",wk);
			continue;
		}
		//if(totWt[i]==0) printf("Erroreeeeeee\n");
		if(mywt[0]==0) printf("Erroreeeeeee\n");
		fflush(stdout);
		//endrun(4321);
	        //wk_p[0]=wk;
		wk/=mywt[0];
#ifdef WINDS
		if(mode==1) //Winds
		{
		 e_shock=0;//edotWIND*All.TimeStep;
		 p_shock=pdotWIND*All.TimeStep;
		 m_shock=mdotWIND*All.TimeStep;
		 met_shock=P[iSN].Metallicity[0]*m_shock;
		}
		else //SNe
#endif
		{
		 e_shock=eSN;
		 p_shock=0;
		 m_shock=mSN;
		 met_shock=zSN;
		}

		//printf("%g %g %g\n",mdotWIND,m_shock,met_shock);
		SphP[j].InternalEnergy += e_shock*wk/P[j].Mass;

#ifdef WINDS
		//printf("%g %g -> %g\n",All.xSNe[iSN],P[j].Pos[0],wk_p[0]);
		double x=P[j].Mass/(P[j].Mass+m_shock*wk);
		//double v=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]+P[j].Vel[2]*P[j].Vel[2]);
		for(k=0;k<3;k++)
		{
		  wk_p[k]/=mywt[1];
		  //wk_p[k+3]/=mywt[k+4];
		  //Momentum
		  P[j].Vel[k]        =x*(P[j].Vel[k]        - p_shock*(wk_p[k])/P[j].Mass);
		  SphP[j].VelPred[k] =x*(SphP[j].VelPred[k] - p_shock*(wk_p[k])/P[j].Mass);
		}
		//double vnew=sqrt(P[j].Vel[0]*P[j].Vel[0]+P[j].Vel[1]*P[j].Vel[1]+P[j].Vel[2]*P[j].Vel[2]);
		//printf("Delta v: %g / Delta p: %g of %g\n",vnew-v,P[j].Mass*(vnew-v),p_shock);
		//fflush(stdout);
		//tot_mom+=p_shock*sqrt(pow(wk_p[0],2)+pow(wk_p[1],2)+pow(wk_p[2],2));
		P[j].Metallicity[0]=x*(P[j].Metallicity[0]+met_shock*wk/P[j].Mass);
		if(P[j].Mass<m_shock*wk)
		{  
			printf("Mass: %d %g %d %g %g %g %g %g\n",iSN,m_shock,numngb_inbox,mywt[0],P[j].Mass,sqrt(r2),wk,wk_p[0]);
			fflush(stdout);
			getchar();
			endrun(1234);
		}
		P[j].Mass+=m_shock*wk;
		//out[i]+=wk;
#endif

              } // for(n = 0; n < numngb; n++)
		//if(iSN==7 && mywt>0) printf("particle %d: %g\n",ThisTask,mywt);		
          } // while(startnode >= 0)
        
        } // while(startnode >= 0)
    //if(iSN==7 && mywt>0)	printf("particle out%d: %g %d\n",iSN,mywt,ThisTask);		
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("My tot: %d %g %d\n",mode,mywt,ThisTask);
    //printf("My (%d) tot momentum for %d is %g of %g\n",ThisTask,iSN, tot_mom,p_shock);
    //fflush(stdout);
    for(k=0;k<NWT;k++)
	outw[k]=mywt[k];
    //getchar();
    myfree(Ngblist);
} // void loop_neighbour
