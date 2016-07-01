#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif

#include "allvars.h"
#include "proto.h"
#include "system/vector.h"


/*! \file ngb.c
 *  \brief neighbour search by means of the tree
 *
 *  This file contains routines for neighbour finding.  We use the
 *  gravity-tree and a range-searching technique to find neighbours.
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * in part by Phil Hopkins (phopkins@caltech.edu) for GIZMO (adding/consolidating 
 * some of the search routines as needed for different fluids)
 */


#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport, mutex_partnodedrift;

#define LOCK_NEXPORT         pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT       pthread_mutex_unlock(&mutex_nexport);
#define LOCK_PARTNODEDRIFT   pthread_mutex_lock(&mutex_partnodedrift);
#define UNLOCK_PARTNODEDRIFT pthread_mutex_unlock(&mutex_partnodedrift);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT
#endif


#ifdef DO_NOT_BRACH_IF

#ifdef __xlC__
#pragma alloca
#define ALLOC_STACK(n) alloca(n)
#elif defined(__GNUC__)
#define ALLOC_STACK(n) alloca(n)
#elif defined(__INTEL_COMPILER)
#define ALLOC_STACK(n) _alloca(n)
#else
#define ALLOC_STACK(n) alloca(n)
#endif


int ngb_filter_pairs(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
		     MyFloat hsml) __attribute__ ((noinline));
int ngb_filter_variables(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
			 MyFloat hsml) __attribute__ ((noinline));

int ngb_filter_pairs(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
		     MyFloat hsml)
{
  int numngb_old = numngb;
  int *comp;
  int no;
  if(!(comp = ALLOC_STACK(numngb_old*sizeof(long long))))
    {
      printf("Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
	     numngb_old*sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  // comp = (int *) mymalloc("NgbFilter", numngb_old * sizeof(int));

  // first compute all the distances                                                                                                                                            

  numngb = 0;
  for(no = 0; no < numngb_old; no++)
    {
      int p = list[no];
      MyDouble dx, dy, dz, d2;
      MyDouble dist = DMAX(PPP[p].Hsml, hsml);

        dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
        dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
        dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
      d2 = dx * dx + dy * dy + dz * dz;
      comp[no] = (d2 < dist * dist);
    }

  // then filter out the distant particles
  if(numngb_old > 0)
    for(no = 0; no < numngb_old; no++)
      {
	if(comp[no])
	  list[numngb++] = list[no];
      }

  //  myfree(comp);

#ifdef NGB_DEBUG
  printf("ngb_treefind_pairs: numngb before/after filter: %d / %d\n", numngb_old, numngb);
#endif
#
  return (int) numngb;
}

int ngb_filter_variables(long long numngb, int list[], t_vector * center, t_vector * box, t_vector * hbox,
			 MyFloat dist)
{
  int numngb_old = numngb;
  int *comp;
  int no;
  if(!(comp = ALLOC_STACK(numngb_old*sizeof(long long))))
    {
      printf("Failed to allocate additional memory for `comp' (%d Mbytes), switch off 'DO_NOT_BRACH_IF'.\n",
	     numngb_old*sizeof(long long));
      endrun(124);
    }
  // mymalloc is nod thread save !!
  //  comp = (int *) mymalloc("NgbFilter", numngb_old * sizeof(int));

  // first compute all the distances                                                                                                                                            

  numngb = 0;
  for(no = 0; no < numngb_old; no++)
    {
      int p = list[no];
      MyDouble dx, dy, dz, d2;

        dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
        dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
        dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - center->d[0], P[p].Pos[1] - center->d[1], P[p].Pos[2] - center->d[2],-1);
      d2 = dx * dx + dy * dy + dz * dz;
      comp[no] = (d2 < dist * dist);
    }

  // then filter out the distant particles
  if(numngb_old > 0)
    for(no = 0; no < numngb_old; no++)
      {
	if(comp[no])
	  list[numngb++] = list[no];
      }

  //  myfree(comp);

#ifdef NGB_DEBUG
  printf("ngb_treefind_vars: numngb before/after filter: %d / %d\n", numngb_old, numngb);
#endif
#
  return numngb;
}


#endif // DO_NOT_BRACH_IF


/*! This routine finds all neighbours `j' that can interact with the
 *  particle `i' in the communication buffer.
 *
 *  Note that an interaction can take place if 
 *  \f$ r_{ij} < h_i \f$  OR if  \f$ r_{ij} < h_j \f$. 
 * 
 *  In the range-search this is taken into account, i.e. it is guaranteed that
 *  all particles are found that fulfil this condition, including the (more
 *  difficult) second part of it. For this purpose, each node knows the
 *  maximum h occuring among the particles it represents.
 */
int ngb_treefind_pairs(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
		       int mode, int *nexport, int *nsend_local)
{
  int no, p, numngb, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;
#ifndef DO_NOT_BRACH_IF
  MyDouble dist, dx, dy, dz;
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  nexport_save = *nexport;

  numngb = 0;

  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  if(P[p].Mass <= 0)
	    continue;

	  if(P[p].Ti_current != ti_Current)
	    drift_particle(p, ti_Current);

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(PPP[p].Hsml, hsml);

	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif
	  Ngblist[numngb++] = p;	/* Note: unlike in previous versions of the code, the buffer 
					   can hold up to all particles */
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(23131);

	      if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= bunchSize)
		    {
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(13003);	/* in this case, the buffer is too small to process even a single particle */
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;
		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (maxPart + maxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_pairs(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    force_drift_node(no, ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len;

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, DMAX(Extnodes[no].hmax, hsml));
#endif
	}
    }


  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_pairs(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}


int ngb_treefind_pairs_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
			       int mode, int *exportflag, int *exportnodecount, int *exportindex,
			       int *ngblist)
{
  int no, p, numngb, task, nexp;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;
#ifndef DO_NOT_BRACH_IF
  MyDouble dist, dx, dy, dz;
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  numngb = 0;

  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

        if(P[p].Type > 0)
            continue;
        
        if(P[p].Mass <= 0)
            continue;

	  if(P[p].Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(PPP[p].Hsml, hsml);

	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif

	  ngblist[numngb++] = p;	/* Note: unlike in previous versions of the code, the buffer 
					   can hold up to all particles */
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
#ifdef DONOTUSENODELIST
	      if(mode == 1)
		{
		  no = Nextnode[no - maxNodes];
		  continue;
		}
#endif
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
		      int exitFlag = 0;
		      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
		      {
			if(Nexport >= bunchSize)
			  {
			    /* out of buffer space. Need to discard work for this particle and interrupt */
			    BufferFullFlag = 1;
			    exitFlag = 1;
			  }
			else
			  {
			    nexp = Nexport;
			    Nexport++;
			  }
		      }
		      UNLOCK_NEXPORT;
		      if(exitFlag)
			return -1;

		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

#ifndef DONOTUSENODELIST
		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
#endif
		}

	      no = Nextnode[no - maxNodes];
	      continue;

	    }

	  current = &Nodes[no];

#ifndef DONOTUSENODELIST
	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }
#endif

	  if(current->Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      force_drift_node(no, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len;

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, DMAX(Extnodes[no].hmax, hsml));
#endif
	}
    }


  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_pairs(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}





/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
int ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			  int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;
#ifndef DO_NOT_BRACH_IF
  MyDouble dist, dx, dy, dz;
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

        if(P[p].Type > 0)
            continue;
        
        if(P[p].Mass <= 0)
            continue;

	  if(P[p].Ti_current != ti_Current)
	    drift_particle(p, ti_Current);

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;

	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= bunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13004);	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++)
			    nsend_local[task] = 0;
			  for(no = 0; no < nexport_save; no++)
			    nsend_local[DataIndexTable[no].Task]++;
			  return -1;
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    force_drift_node(no, ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }


  //printf("%d\n", numngb);

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}


/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				  int mode, int *exportflag, int *exportnodecount, int *exportindex,
				  int *ngblist)
{
  int numngb, no, nexp, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;

#ifndef DO_NOT_BRACH_IF
  MyDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

        if(P[p].Type > 0)
            continue;
        
        if(P[p].Mass <= 0)
            continue;

	  if(P[p].Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      drift_particle(p, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif
	  ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      exportflag[task] = target;
		      exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(exportnodecount[task] == NODELISTLENGTH)
		    {
		      int exitFlag = 0;
		      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
		      {
			if(Nexport >= bunchSize)
			  {
			    /* out of buffer space. Need to discard work for this particle and interrupt */
			    BufferFullFlag = 1;
			    exitFlag = 1;
			  }
			else
			  {
			    nexp = Nexport;
			    Nexport++;
			  }
		      }
		      UNLOCK_NEXPORT;
		      if(exitFlag)
			return -1;

		      exportnodecount[task] = 0;
		      exportindex[task] = nexp;
		      DataIndexTable[nexp].Task = task;
		      DataIndexTable[nexp].Index = target;
		      DataIndexTable[nexp].IndexGet = nexp;
		    }

		  DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;

		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(current->Ti_current != ti_Current)
	    {
	      LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
	      force_drift_node(no, ti_Current);
	      UNLOCK_PARTNODEDRIFT;
	    }

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}





/*! This function constructs the neighbour tree. To this end, we actually need
 *  to construct the gravitational tree, because we use it now for the
 *  neighbour search.
 */
void ngb_treebuild(void)
{
  if(ThisTask == 0)
    printf("Begin Ngb-tree construction.\n");

  CPU_Step[CPU_MISC] += measure_time();

  force_treebuild(NumPart, NULL);

  CPU_Step[CPU_TREEBUILD] += measure_time();

  if(ThisTask == 0)
    printf("Ngb-Tree contruction finished \n");
}



int ngb_treefind_fof_primary(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  MyDouble dx, dy, dz, dist, r2;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;

#ifndef DO_NOT_BRACH_IF
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(!((1 << P[p].Type) & (MyFLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES)))
	    continue;

	  if(mode == 0)
	    continue;

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif
	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(mode == 0)
		{
		  if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= bunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++)
			    nsend_local[task] = 0;
			  for(no = 0; no < nexport_save; no++)
			    nsend_local[DataIndexTable[no].Task]++;
			  return -1;
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      if(mode == -1)
		{
		  *nexport = 1;
		}

	      no = Nextnode[no - maxNodes];
	      continue;

	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

	  if(mode == 0)
	    {
	      if(!(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL)))	/* we have a node with only local particles, can skip branch */
		{
		  no = current->u.d.sibling;
		  continue;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
	    continue;

	  if((current->u.d.bitflags & ((1 << BITFLAG_TOPLEVEL) + (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))) == 0)	/* only use fully local nodes */
	    {
	      /* test whether the node is contained within the sphere */
	      dist = hsml - FACT2 * current->len;
	      if(dist > 0)
		if(r2 < dist * dist)
		  {
		    if(current->u.d.bitflags & (1 << BITFLAG_INSIDE_LINKINGLENGTH))	/* already flagged */
		      {
			/* sufficient to return only one particle inside this cell */

			p = current->u.d.nextnode;
			while(p >= 0)
			  {
			    if(p < maxPart)
			      {
				if(((1 << P[p].Type) & (MyFLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES)))
				  {
#ifndef DO_NOT_BRACH_IF
				    dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
				    dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
				    dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
				    if(dx * dx + dy * dy + dz * dz > hsml * hsml)
				      break;
#endif
				    Ngblist[numngb++] = p;
				    break;
				  }
				p = Nextnode[p];
			      }
			    else if(p >= maxPart + maxNodes)
			      p = Nextnode[p - maxNodes];
			    else
			      p = Nodes[p].u.d.nextnode;
			  }
			continue;
		      }
		    else
		      {
			/* flag it now */
			current->u.d.bitflags |= (1 << BITFLAG_INSIDE_LINKINGLENGTH);
		      }
		  }
	    }

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}



/* find all particles of type FLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES in kernel length in order to find nearest one */
int ngb_treefind_fof_nearest(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
#ifndef DO_NOT_BRACH_IF
  MyDouble dx, dy, dz, dist;
#ifdef PERIODIC
  MyDouble xtmp;
#endif
#else
  t_vector box, hbox, vcenter;
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif

#define FACT2 0.86602540

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(!((1 << P[p].Type) & (MyFLAG_NOT_IN_PUBLIC_CODE_PRIMARY_LINK_TYPES)))
	    continue;

#ifndef DO_NOT_BRACH_IF
	  dist = hsml;
	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;
#endif

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(123192);

	      if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		{
		  Exportflag[task] = target;
		  Exportnodecount[task] = NODELISTLENGTH;
		}

	      if(Exportnodecount[task] == NODELISTLENGTH)
		{
		  if(*nexport >= bunchSize)
		    {
		      *nexport = nexport_save;
		      if(nexport_save == 0)
			endrun(13005);	/* in this case, the buffer is too small to process even a single particle */
		      for(task = 0; task < NTask; task++)
			nsend_local[task] = 0;
		      for(no = 0; no < nexport_save; no++)
			nsend_local[DataIndexTable[no].Task]++;
		      return -1;
		    }
		  Exportnodecount[task] = 0;
		  Exportindex[task] = *nexport;
		  DataIndexTable[*nexport].Task = task;
		  DataIndexTable[*nexport].Index = target;
		  DataIndexTable[*nexport].IndexGet = *nexport;
		  *nexport = *nexport + 1;
		  nsend_local[task]++;
		}

	      DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		DomainNodeIndex[no - (maxPart + maxNodes)];

	      if(Exportnodecount[task] < NODELISTLENGTH)
		DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
		  return numngb;
#else
		  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
		}
	    }

#ifndef DO_NOT_BRACH_IF
	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;;
	  dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if((dx * dx + dy * dy + dz * dz) > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
	  no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
	}
    }

  *startnode = -1;
#ifndef DO_NOT_BRACH_IF
  return numngb;
#else
  return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml);
#endif
}



















#ifdef ADAPTIVE_GRAVSOFT_FORALL
/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
int ags_ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
                                  int mode, int *exportflag, int *exportnodecount, int *exportindex,
                                  int *ngblist, int type_of_searching_particle)
{
    int numngb, no, nexp, p, task;
    struct NODE *current;
    // cache some global vars locally for improved compiler alias analysis
    int maxPart = All.MaxPart;
    int maxNodes = MaxNodes;
    int bunchSize = All.BunchSize;
    integertime ti_Current = All.Ti_Current;
    
#ifndef DO_NOT_BRACH_IF
    MyDouble dx, dy, dz, dist;
#ifdef PERIODIC
    MyDouble xtmp;
#endif
#else
    t_vector box, hbox, vcenter;
    INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
    INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
    SCALE_VECTOR3(0.5, &box, &hbox);
#endif
    

    numngb = 0;
    no = *startnode;
    
    while(no >= 0)
    {
        if(no < maxPart)		/* single particle */
        {
            p = no;
            no = Nextnode[no];
            
            /* call the master routine which decides if particles "talk to" each other in-kernel */
            if(ags_gravity_kernel_shared_check(type_of_searching_particle , P[p].Type) == 0)
                continue;
            
            if(P[p].Mass <= 0)
                continue;
            
            if(P[p].Ti_current != ti_Current)
            {
                LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                drift_particle(p, ti_Current);
                UNLOCK_PARTNODEDRIFT;
            }
            
#ifndef DO_NOT_BRACH_IF
            dist = hsml;
            dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
            if(dz > dist)
                continue;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
#endif
            ngblist[numngb++] = p;
        }
        else
        {
            if(no >= maxPart + maxNodes)	/* pseudo particle */
            {
                if(mode == 1)
                    endrun(12312);
                
                if(target >= 0)	/* if no target is given, export will not occur */
                {
                    if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
                    {
                        exportflag[task] = target;
                        exportnodecount[task] = NODELISTLENGTH;
                    }
                    
                    if(exportnodecount[task] == NODELISTLENGTH)
                    {
                        int exitFlag = 0;
                        LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
                        {
                            if(Nexport >= bunchSize)
                            {
                                /* out of buffer space. Need to discard work for this particle and interrupt */
                                BufferFullFlag = 1;
                                exitFlag = 1;
                            }
                            else
                            {
                                nexp = Nexport;
                                Nexport++;
                            }
                        }
                        UNLOCK_NEXPORT;
                        if(exitFlag)
                            return -1;
                        
                        exportnodecount[task] = 0;
                        exportindex[task] = nexp;
                        DataIndexTable[nexp].Task = task;
                        DataIndexTable[nexp].Index = target;
                        DataIndexTable[nexp].IndexGet = nexp;
                    }
                    
                    DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] =
                    DomainNodeIndex[no - (maxPart + maxNodes)];
                    
                    if(exportnodecount[task] < NODELISTLENGTH)
                        DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
                    
                }
                
                no = Nextnode[no - maxNodes];
                continue;
            }
            
            current = &Nodes[no];
            
            if(mode == 1)
            {
                if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                {
                    *startnode = -1;
#ifndef DO_NOT_BRACH_IF
                    return numngb;
#else
                    return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
                }
            }
            
            if(current->Ti_current != ti_Current)
            {
                LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                force_drift_node(no, ti_Current);
                UNLOCK_PARTNODEDRIFT;
            }
            
            if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
            {
                if(current->u.d.mass)	/* open cell */
                {
                    no = current->u.d.nextnode;
                    continue;
                }
            }
            
#ifndef DO_NOT_BRACH_IF
            no = current->u.d.sibling;	/* in case the node can be discarded */
            
            dist = hsml + 0.5 * current->len;
            dx = NGB_PERIODIC_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dx > dist)
                continue;
            dy = NGB_PERIODIC_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dy > dist)
                continue;
            dz = NGB_PERIODIC_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
            if(dz > dist)
                continue;
            /* now test against the minimal sphere enclosing everything */
            dist += FACT1 * current->len;
            if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;
            
            no = current->u.d.nextnode;	/* ok, we need to open the node */
#else
            no = ngb_check_node(current, &vcenter, &box, &hbox, hsml);
#endif
        }
    }
    
    *startnode = -1;
#ifndef DO_NOT_BRACH_IF
    return numngb;
#else
    return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml);
#endif
}

#endif // ADAPTIVE_GRAVSOFT_FORALL









