#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif


/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). The code has been modified
 * substantially (condensed, new feedback routines added,
 * some optimizatins, and new variable/memory conventions added)
 * by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

/*! auxiliary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float shortrange_table[NTAB], shortrange_table_potential[NTAB];
/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;

static int tree_allocated_flag = 0;


#ifdef OMP_NUM_THREADS
extern pthread_mutex_t mutex_nexport, mutex_partnodedrift, mutex_workcount;

#define LOCK_NEXPORT         pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT       pthread_mutex_unlock(&mutex_nexport);
#define LOCK_PARTNODEDRIFT   pthread_mutex_lock(&mutex_partnodedrift);
#define UNLOCK_PARTNODEDRIFT pthread_mutex_unlock(&mutex_partnodedrift);

/*! The cost computation for the tree-gravity (required for the domain
 decomposition) is not exactly thread-safe if THREAD_SAFE_COSTS is not defined.
 However using locks for an exactly thread-safe cost computiation results in a
 significant (~25%) performance penalty in the tree-walk while having only an
 extremely small effect on the obtained costs. The domain decomposition should
 thus not be significantly changed if THREAD_SAFE_COSTS is not used.*/
#ifdef THREAD_SAFE_COSTS
#define LOCK_WORKCOUNT       pthread_mutex_lock(&mutex_workcount);
#define UNLOCK_WORKCOUNT     pthread_mutex_unlock(&mutex_workcount);
#else
#define LOCK_WORKCOUNT
#define UNLOCK_WORKCOUNT
#endif

#else

#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#define LOCK_PARTNODEDRIFT
#define UNLOCK_PARTNODEDRIFT
#define LOCK_WORKCOUNT
#define UNLOCK_WORKCOUNT

#endif



#ifdef PERIODIC
/*! Size of 3D lock-up table for Ewald correction force */
#define EN  64
/*! 3D lock-up table for Ewald correction to force and potential. Only one
 *  octant is stored, the rest constructed by using the symmetry
 */
static MyFloat fcorrx[EN + 1][EN + 1][EN + 1];
static MyFloat fcorry[EN + 1][EN + 1][EN + 1];
static MyFloat fcorrz[EN + 1][EN + 1][EN + 1];
static MyFloat potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
#endif



/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart, struct unbind_data *mp)
{
    int flag;
    
    
    do
    {
        Numnodestree = force_treebuild_single(npart, mp);
        
        MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if(flag == -1)
        {
            force_treefree();
            
            if(ThisTask == 0)
                printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);
            
            All.TreeAllocFactor *= 1.15;
            
            if(ThisTask == 0)
                printf("new value=%g\n", All.TreeAllocFactor);
            
            force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
        }
    }
    while(flag == -1);
    
    force_flag_localnodes();
    
    force_exchange_pseudodata();
    
    force_treeupdate_pseudos(All.MaxPart);
    
    TimeOfLastTreeConstruction = All.Time;
    
    return Numnodestree;
}



/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart, struct unbind_data *mp)
{
    int i, j, k, subnode = 0, shift, parent, numnodes, rep;
    int nfree, th, nn, no;
    struct NODE *nfreep;
    MyFloat lenhalf;
    peanokey key, morton, th_key, *morton_list;
    
    
    /* create an empty root node  */
    nfree = All.MaxPart;		/* index of first free node */
    nfreep = &Nodes[nfree];	/* select first node */
    
    nfreep->len = DomainLen;
    for(j = 0; j < 3; j++)
        nfreep->center[j] = DomainCenter[j];
    for(j = 0; j < 8; j++)
        nfreep->u.suns[j] = -1;
    
    
    numnodes = 1;
    nfreep++;
    nfree++;
    
    /* create a set of empty nodes corresponding to the top-level domain
     * grid. We need to generate these nodes first to make sure that we have a
     * complete top-level tree which allows the easy insertion of the
     * pseudo-particles at the right place
     */
    
    force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);
    
    /* if a high-resolution region in a global tree is used, we need to generate
     * an additional set empty nodes to make sure that we have a complete
     * top-level tree for the high-resolution inset
     */
    
    nfreep = &Nodes[nfree];
    parent = -1;			/* note: will not be used below before it is changed */
    
    morton_list = (peanokey *) mymalloc("morton_list", NumPart * sizeof(peanokey));
    
    /* now we insert all particles */
    for(k = 0; k < npart; k++)
    {
        if(mp)
            i = mp[k].index;
        else
            i = k;
        
#ifdef NEUTRINOS
        if(P[i].Type == 2)
            continue;
#endif
        
        rep = 0;
        
        key = peano_and_morton_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
                                   (int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
                                   (int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION,
                                   &morton);
        morton_list[i] = morton;
        
        shift = 3 * (BITS_PER_DIMENSION - 1);
        
        no = 0;
        while(TopNodes[no].Daughter >= 0)
        {
            no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);
            shift -= 3;
            rep++;
        }
        
        no = TopNodes[no].Leaf;
        th = DomainNodeIndex[no];
        
        while(1)
        {
            if(th >= All.MaxPart)	/* we are dealing with an internal node */
            {
                if(shift >= 0)
                {
                    subnode = ((morton >> shift) & 7);
                }
                else
                {
                    subnode = 0;
                    if(P[i].Pos[0] > Nodes[th].center[0])
                        subnode += 1;
                    if(P[i].Pos[1] > Nodes[th].center[1])
                        subnode += 2;
                    if(P[i].Pos[2] > Nodes[th].center[2])
                        subnode += 4;
                }
                
#ifndef NOTREERND
                if(Nodes[th].len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[P[i].Type])
                {
                    /* seems like we're dealing with particles at identical (or extremely close)
                     * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening
                     * length-scale anyway.
                     */
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
                    subnode = (int) (8.0 * get_random_number((P[i].ID + rep) % (RNDTABLE + (rep & 3))));
#else
                    subnode = (int) (8.0 * get_random_number(P[i].ID));
#endif
                    
                    if(subnode >= 8)
                        subnode = 7;
                }
#endif
                
                nn = Nodes[th].u.suns[subnode];
                
                shift -= 3;
                
                if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
                {
                    parent = th;
                    th = nn;
                    rep++;
                }
                else
                {
                    /* here we have found an empty slot where we can attach
                     * the new particle as a leaf.
                     */
                    Nodes[th].u.suns[subnode] = i;
                    break;	/* done for this particle */
                }
            }
            else
            {
                /* We try to insert into a leaf with a single particle.  Need
                 * to generate a new internal node at this point.
                 */
                Nodes[parent].u.suns[subnode] = nfree;
                
                nfreep->len = 0.5 * Nodes[parent].len;
                lenhalf = 0.25 * Nodes[parent].len;
                
                if(subnode & 1)
                    nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
                else
                    nfreep->center[0] = Nodes[parent].center[0] - lenhalf;
                
                if(subnode & 2)
                    nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
                else
                    nfreep->center[1] = Nodes[parent].center[1] - lenhalf;
                
                if(subnode & 4)
                    nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
                else
                    nfreep->center[2] = Nodes[parent].center[2] - lenhalf;
                
                nfreep->u.suns[0] = -1;
                nfreep->u.suns[1] = -1;
                nfreep->u.suns[2] = -1;
                nfreep->u.suns[3] = -1;
                nfreep->u.suns[4] = -1;
                nfreep->u.suns[5] = -1;
                nfreep->u.suns[6] = -1;
                nfreep->u.suns[7] = -1;
                
                if(shift >= 0)
                {
                    th_key = morton_list[th];
                    subnode = ((th_key >> shift) & 7);
                }
                else
                {
                    subnode = 0;
                    if(P[th].Pos[0] > nfreep->center[0])
                        subnode += 1;
                    if(P[th].Pos[1] > nfreep->center[1])
                        subnode += 2;
                    if(P[th].Pos[2] > nfreep->center[2])
                        subnode += 4;
                }
                
#ifndef NOTREERND
                if(nfreep->len < EPSILON_FOR_TREERND_SUBNODE_SPLITTING * All.ForceSoftening[P[th].Type])
                {
                    /* seems like we're dealing with particles at identical (or extremely close)
                     * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but this will only happen well below gravitational softening
                     * length-scale anyway.
                     */
#ifdef USE_PREGENERATED_RANDOM_NUMBER_TABLE
                    subnode = (int) (8.0 * get_random_number((P[th].ID + rep) % (RNDTABLE + (rep & 3))));
#else
                    subnode = (int) (8.0 * get_random_number(P[th].ID));
#endif
                    
                    if(subnode >= 8)
                        subnode = 7;
                }
#endif
                nfreep->u.suns[subnode] = th;
                
                th = nfree;	/* resume trying to insert the new particle at
                             * the newly created internal node
                             */
                
                numnodes++;
                nfree++;
                nfreep++;
                
                if((numnodes) >= MaxNodes)
                {
                    printf("task %d: maximum number %d of tree-nodes reached for particle %d.\n", ThisTask,
                           MaxNodes, i);
                    
                    if(All.TreeAllocFactor > 5.0)
                    {
                        printf
                        ("task %d: looks like a serious problem for particle %d, stopping with particle dump.\n",
                         ThisTask, i);
                        dump_particles();
                        endrun(1);
                    }
                    else
                    {
                        myfree(morton_list);
                        return -1;
                    }
                }
            }
        }
    }
    
    myfree(morton_list);
    
    
    /* insert the pseudo particles that represent the mass distribution of other domains */
    force_insert_pseudo_particles();
    
    
    /* now compute the multipole moments recursively */
    last = -1;
    
    force_update_node_recursive(All.MaxPart, -1, -1);
    
    if(last >= All.MaxPart)
    {
        if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
            Nextnode[last - MaxNodes] = -1;
        else
            Nodes[last].u.d.nextnode = -1;
    }
    else
        Nextnode[last] = -1;
    
    return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
                              int *nextfree)
{
    int i, j, k, n, sub, count;
    MyFloat lenhalf;
    
    if(TopNodes[topnode].Daughter >= 0)
    {
        for(i = 0; i < 2; i++)
            for(j = 0; j < 2; j++)
                for(k = 0; k < 2; k++)
                {
                    sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);
                    
                    count = i + 2 * j + 4 * k;
                    
                    Nodes[no].u.suns[count] = *nextfree;
                    
                    lenhalf = 0.25 * Nodes[no].len;
                    Nodes[*nextfree].len = 0.5 * Nodes[no].len;
                    Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
                    Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
                    Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;
                    
                    for(n = 0; n < 8; n++)
                        Nodes[*nextfree].u.suns[n] = -1;
                    
                    if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                        DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;
                    
                    *nextfree = *nextfree + 1;
                    *nodecount = *nodecount + 1;
                    
                    if((*nodecount) >= MaxNodes || (*nodecount) >= MaxTopNodes)
                    {
                        printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
                               "MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
                               ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
                        printf("in create empty nodes\n");
                        dump_particles();
                        endrun(11);
                    }
                    
                    force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
                                             bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
                }
    }
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
    int i, index;
    
    for(i = 0; i < NTopleaves; i++)
    {
        index = DomainNodeIndex[i];
        
        if(DomainTask[i] != ThisTask)
            Nodes[index].u.suns[0] = All.MaxPart + MaxNodes + i;
    }
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 */
void force_update_node_recursive(int no, int sib, int father)
{
    int j, jj, k, p, pp, nextsib, suns[8], count_particles, multiple_flag;
    MyFloat hmax, vmax, v;
    MyFloat divVmax, divVel;
    MyFloat s[3], vs[3], mass;
    struct particle_data *pa;
    
#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS], sigma_eff=0;
    for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
    
    MyFloat maxsoft;
    
    if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
        for(j = 0; j < 8; j++)
            suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will overwrite one element (union!) */
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
                    Nextnode[last - MaxNodes] = no;
                else
                    Nodes[last].u.d.nextnode = no;
            }
            else
                Nextnode[last] = no;
        }
        
        last = no;
        
#ifdef RT_USE_GRAVTREE
        for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
        mass = 0;
        s[0] = 0;
        s[1] = 0;
        s[2] = 0;
        vs[0] = 0;
        vs[1] = 0;
        vs[2] = 0;
        hmax = 0;
        vmax = 0;
        divVmax = 0;
        count_particles = 0;
        maxsoft = 0;
        
        for(j = 0; j < 8; j++)
        {
            if((p = suns[j]) >= 0)
            {
                /* check if we have a sibling on the same level */
                for(jj = j + 1; jj < 8; jj++)
                    if((pp = suns[jj]) >= 0)
                        break;
                
                if(jj < 8)	/* yes, we do */
                    nextsib = pp;
                else
                    nextsib = sib;
                
                force_update_node_recursive(p, nextsib, no);
                
                if(p >= All.MaxPart)	/* an internal node or pseudo particle */
                {
                    if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
                    {
                        /* nothing to be done here because the mass of the
                         * pseudo-particle is still zero. This will be changed
                         * later.
                         */
                    }
                    else
                    {
                        mass += (Nodes[p].u.d.mass);
                        s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
                        s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
                        s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
                        vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
                        vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
                        vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);
#ifdef RT_USE_GRAVTREE
                        for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += (Nodes[p].stellar_lum[k]);}
#endif
                        if(Nodes[p].u.d.mass > 0)
                        {
                            if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
                                count_particles += 2;
                            else
                                count_particles++;
                        }
                        
                        if(Extnodes[p].hmax > hmax)
                            hmax = Extnodes[p].hmax;
                        
                        if(Extnodes[p].vmax > vmax)
                            vmax = Extnodes[p].vmax;
                        if(Extnodes[p].divVmax > divVmax)
                            divVmax = Extnodes[p].divVmax;
                        
                        /* update of the maximum gravitational softening in the node */
                        if(Nodes[p].maxsoft > maxsoft)
                            maxsoft = Nodes[p].maxsoft;
                        
                    }
                }
                else		/* a particle */
                {
                    count_particles++;
                    
                    pa = &P[p];
                    
                    mass += (pa->Mass);
                    s[0] += (pa->Mass * pa->Pos[0]);
                    s[1] += (pa->Mass * pa->Pos[1]);
                    s[2] += (pa->Mass * pa->Pos[2]);
                    vs[0] += (pa->Mass * pa->Vel[0]);
                    vs[1] += (pa->Mass * pa->Vel[1]);
                    vs[2] += (pa->Mass * pa->Vel[2]);
                    
#ifdef RT_USE_GRAVTREE
                    double lum[N_RT_FREQ_BINS];
                    int active_check = rt_get_source_luminosity(p,sigma_eff,lum);
                    if(active_check)
                    {
                        double l_sum = 0; for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += lum[k]; l_sum += lum[k];}
                    }
#endif
                    
                    
                    
                    
                    
                    if(pa->Type == 0)
                    {
                        if(PPP[p].Hsml > hmax)
                            hmax = PPP[p].Hsml;
                        
                        divVel = P[p].Particle_DivVel;
                        if(divVel > divVmax)
                            divVmax = divVel;
                    }
                    
                    for(k = 0; k < 3; k++)
                        if((v = fabs(pa->Vel[k])) > vmax)
                            vmax = v;
                    
                    /* update of the maximum gravitational softening  */
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                    if(PPP[p].AGS_Hsml > maxsoft)
                        maxsoft = PPP[p].AGS_Hsml;
#else
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
                    if(All.ForceSoftening[pa->Type] > maxsoft)
                        maxsoft = All.ForceSoftening[pa->Type];
#else
                    if(pa->Type == 0)
                    {
                        if(PPP[p].Hsml > maxsoft)
                            maxsoft = PPP[p].Hsml;
                    }
                    else
                    {
                        if(All.ForceSoftening[pa->Type] > maxsoft)
                            maxsoft = All.ForceSoftening[pa->Type];
                    }
#endif
#endif
                }
            }
        }
        
        
        if(mass)
        {
            s[0] /= mass;
            s[1] /= mass;
            s[2] /= mass;
            vs[0] /= mass;
            vs[1] /= mass;
            vs[2] /= mass;
        }
        else
        {
            s[0] = Nodes[no].center[0];
            s[1] = Nodes[no].center[1];
            s[2] = Nodes[no].center[2];
            vs[0] = 0;
            vs[1] = 0;
            vs[2] = 0;
        }
        
        
        
        Nodes[no].Ti_current = All.Ti_Current;
        Nodes[no].u.d.mass = mass;
        Nodes[no].u.d.s[0] = s[0];
        Nodes[no].u.d.s[1] = s[1];
        Nodes[no].u.d.s[2] = s[2];
        Nodes[no].GravCost = 0;
#ifdef RT_USE_GRAVTREE
        for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = stellar_lum[k];}
#endif
        
        Extnodes[no].Ti_lastkicked = All.Ti_Current;
        Extnodes[no].Flag = GlobFlag;
        Extnodes[no].vs[0] = vs[0];
        Extnodes[no].vs[1] = vs[1];
        Extnodes[no].vs[2] = vs[2];
        Extnodes[no].hmax = hmax;
        Extnodes[no].vmax = vmax;
        Extnodes[no].divVmax = divVmax;
        Extnodes[no].dp[0] = 0;
        Extnodes[no].dp[1] = 0;
        Extnodes[no].dp[2] = 0;
        
        if(count_particles > 1)	/* this flags that the node represents more than one particle */
            multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
        else
            multiple_flag = 0;
        
        Nodes[no].u.d.bitflags = multiple_flag;
        Nodes[no].maxsoft = maxsoft;
        Nodes[no].u.d.sibling = sib;
        Nodes[no].u.d.father = father;
    }
    else				/* single particle or pseudo particle */
    {
        if(last >= 0)
        {
            if(last >= All.MaxPart)
            {
                if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
                    Nextnode[last - MaxNodes] = no;
                else
                    Nodes[last].u.d.nextnode = no;
            }
            else
                Nextnode[last] = no;
        }
        
        last = no;
        
        if(no < All.MaxPart)	/* only set it for single particles */
            Father[no] = father;
    }
}




/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
    int i, no, m, ta, recvTask;
    int *recvcounts, *recvoffset;
    struct DomainNODE
    {
        MyFloat s[3];
        MyFloat vs[3];
        MyFloat mass;
        MyFloat hmax;
        MyFloat vmax;
        MyFloat divVmax;        
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        MyFloat maxsoft;
#endif
#ifdef RT_USE_GRAVTREE
        MyFloat stellar_lum[N_RT_FREQ_BINS];
#endif
        unsigned int bitflags;
#ifdef PAD_STRUCTURES
#ifndef DOUBLEPRECISION
        int pad[5];
#else
#if (DOUBLEPRECISION+0) == 2
        /* mixed precision */
        int pad[5];
#else
        int pad[3];
#endif
#endif				/* DOUBLEPRECISION  */
#endif
    }
    *DomainMoment;
    
    
    DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
        for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
            i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
        {
            no = DomainNodeIndex[i];
            
            /* read out the multipole moments from the local base cells */
            DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
            DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
            DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
            DomainMoment[i].vs[0] = Extnodes[no].vs[0];
            DomainMoment[i].vs[1] = Extnodes[no].vs[1];
            DomainMoment[i].vs[2] = Extnodes[no].vs[2];
            DomainMoment[i].mass = Nodes[no].u.d.mass;
            DomainMoment[i].hmax = Extnodes[no].hmax;
            DomainMoment[i].vmax = Extnodes[no].vmax;
            DomainMoment[i].divVmax = Extnodes[no].divVmax;
            DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
            DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
#ifdef RT_USE_GRAVTREE
            int k; for(k=0;k<N_RT_FREQ_BINS;k++) {DomainMoment[i].stellar_lum[k] = Nodes[no].stellar_lum[k];}
#endif
        }
    
    /* share the pseudo-particle data accross CPUs */
    recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
    recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
    {
        for(recvTask = 0; recvTask < NTask; recvTask++)
        {
            recvcounts[recvTask] =
            (DomainEndList[recvTask * MULTIPLEDOMAINS + m] - DomainStartList[recvTask * MULTIPLEDOMAINS + m] +
             1) * sizeof(struct DomainNODE);
            recvoffset[recvTask] = DomainStartList[recvTask * MULTIPLEDOMAINS + m] * sizeof(struct DomainNODE);
        }
#ifdef USE_MPI_IN_PLACE
        MPI_Allgatherv(MPI_IN_PLACE, recvcounts[ThisTask],
                       MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
#else
        MPI_Allgatherv(&DomainMoment[DomainStartList[ThisTask * MULTIPLEDOMAINS + m]], recvcounts[ThisTask],
                       MPI_BYTE, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);
#endif
    }
    
    myfree(recvoffset);
    myfree(recvcounts);
    
    
    for(ta = 0; ta < NTask; ta++)
        if(ta != ThisTask)
            for(m = 0; m < MULTIPLEDOMAINS; m++)
                for(i = DomainStartList[ta * MULTIPLEDOMAINS + m]; i <= DomainEndList[ta * MULTIPLEDOMAINS + m]; i++)
                {
                    no = DomainNodeIndex[i];
                    
                    Nodes[no].u.d.s[0] = DomainMoment[i].s[0];
                    Nodes[no].u.d.s[1] = DomainMoment[i].s[1];
                    Nodes[no].u.d.s[2] = DomainMoment[i].s[2];
                    Extnodes[no].vs[0] = DomainMoment[i].vs[0];
                    Extnodes[no].vs[1] = DomainMoment[i].vs[1];
                    Extnodes[no].vs[2] = DomainMoment[i].vs[2];
                    Nodes[no].u.d.mass = DomainMoment[i].mass;
                    Extnodes[no].hmax = DomainMoment[i].hmax;
                    Extnodes[no].vmax = DomainMoment[i].vmax;
                    Extnodes[no].divVmax = DomainMoment[i].divVmax;
                    Nodes[no].u.d.bitflags =
                    (Nodes[no].u.d.bitflags & (~BITFLAG_MASK)) | (DomainMoment[i].bitflags & BITFLAG_MASK);
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                    Nodes[no].maxsoft = DomainMoment[i].maxsoft;
#endif
#ifdef RT_USE_GRAVTREE
                    int k; for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = DomainMoment[i].stellar_lum[k];}
#endif
                }
    
    myfree(DomainMoment);
}



/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
    int j, p, count_particles, multiple_flag;
    MyFloat hmax, vmax;
    MyFloat divVmax;
    MyFloat s[3], vs[3], mass;
    
#ifdef RT_USE_GRAVTREE
    MyFloat stellar_lum[N_RT_FREQ_BINS];
#endif
    
    MyFloat maxsoft;
    
#ifdef RT_USE_GRAVTREE
    for(j=0;j<N_RT_FREQ_BINS;j++) {stellar_lum[j]=0;}
#endif
    mass = 0;
    s[0] = 0;
    s[1] = 0;
    s[2] = 0;
    vs[0] = 0;
    vs[1] = 0;
    vs[2] = 0;
    hmax = 0;
    vmax = 0;
    divVmax = 0;
    count_particles = 0;
    maxsoft = 0;
    
    p = Nodes[no].u.d.nextnode;
    
    for(j = 0; j < 8; j++)	/* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
    {
        if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)	/* internal node */
        {
            if(Nodes[p].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
                force_treeupdate_pseudos(p);
            
            mass += (Nodes[p].u.d.mass);
            s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
            s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
            s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
#ifdef RT_USE_GRAVTREE
            int k; for(k=0;k<N_RT_FREQ_BINS;k++) {stellar_lum[k] += (Nodes[p].stellar_lum[k]);}
#endif
            vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
            vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
            vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);
            
            if(Extnodes[p].hmax > hmax)
                hmax = Extnodes[p].hmax;
            if(Extnodes[p].vmax > vmax)
                vmax = Extnodes[p].vmax;
            if(Extnodes[p].divVmax > divVmax)
                divVmax = Extnodes[p].divVmax;
            
            if(Nodes[p].u.d.mass > 0)
            {
                if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
                    count_particles += 2;
                else
                    count_particles++;
            }
            
            if(Nodes[p].maxsoft > maxsoft)
                maxsoft = Nodes[p].maxsoft;
        }
        else
            endrun(6767);		/* may not happen */
        
        p = Nodes[p].u.d.sibling;
    }
    
    if(mass)
    {
        s[0] /= mass;
        s[1] /= mass;
        s[2] /= mass;
        vs[0] /= mass;
        vs[1] /= mass;
        vs[2] /= mass;
    }
    else
    {
        s[0] = Nodes[no].center[0];
        s[1] = Nodes[no].center[1];
        s[2] = Nodes[no].center[2];
        vs[0] = 0;
        vs[1] = 0;
        vs[2] = 0;
    }
    
    
    
    Nodes[no].u.d.s[0] = s[0];
    Nodes[no].u.d.s[1] = s[1];
    Nodes[no].u.d.s[2] = s[2];
    Extnodes[no].vs[0] = vs[0];
    Extnodes[no].vs[1] = vs[1];
    Extnodes[no].vs[2] = vs[2];
    Nodes[no].u.d.mass = mass;
#ifdef RT_USE_GRAVTREE
    int k; for(k=0;k<N_RT_FREQ_BINS;k++) {Nodes[no].stellar_lum[k] = stellar_lum[k];}
#endif
    
    Extnodes[no].hmax = hmax;
    Extnodes[no].vmax = vmax;
    Extnodes[no].divVmax = divVmax;
    Extnodes[no].Flag = GlobFlag;
    
    
    if(count_particles > 1)
        multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
    else
        multiple_flag = 0;
    
    Nodes[no].u.d.bitflags &= (~BITFLAG_MASK);	/* this clears the bits */
    Nodes[no].u.d.bitflags |= multiple_flag;
    Nodes[no].maxsoft = maxsoft;
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
    int no, i, m;
    
    /* mark all top-level nodes */
    
    for(i = 0; i < NTopleaves; i++)
    {
        no = DomainNodeIndex[i];
        
        while(no >= 0)
        {
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
                break;
            
            Nodes[no].u.d.bitflags |= (1 << BITFLAG_TOPLEVEL);
            
            no = Nodes[no].u.d.father;
        }
        
        /* mark also internal top level nodes */
        
        no = DomainNodeIndex[i];
        no = Nodes[no].u.d.father;
        
        while(no >= 0)
        {
            if(Nodes[no].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
                break;
            
            Nodes[no].u.d.bitflags |= (1 << BITFLAG_INTERNAL_TOPLEVEL);
            
            no = Nodes[no].u.d.father;
        }
    }
    
    /* mark top-level nodes that contain local particles */
    
    for(m = 0; m < MULTIPLEDOMAINS; m++)
        for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
            i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
        {
            no = DomainNodeIndex[i];
            
            if(DomainTask[i] != ThisTask)
                endrun(131231231);
            
            while(no >= 0)
            {
                if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))
                    break;
                
                Nodes[no].u.d.bitflags |= (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS);
                
                no = Nodes[no].u.d.father;
            }
        }
}


/*! When a new additional star particle is created, we can put it into the
 *  tree at the position of the spawning gas particle. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed.
 */
void force_add_star_to_tree(int igas, int istar)
{
    int no;
    no = Nextnode[igas];
    Nextnode[igas] = istar;
    Nextnode[istar] = no;
    Father[istar] = Father[igas];
}



/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
/*! The modern version of this routine handles both the PM-grid and non-PM
 *  cases, unlike the previous version (which used two, redundant, algorithms)
 */
/*! In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access panelty (which reduces cache performance) incurred by the
 *  table.
 */
int force_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)
{
    struct NODE *nop = 0;
    int no, nodesinlist, ptype, ninteractions, nexp, task, listindex = 0;
    double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
    double pos_x, pos_y, pos_z, aold;
    MyLongDouble acc_x, acc_y, acc_z;
    // cache some global vars in local vars to help compiler with alias analysis
    int maxPart = All.MaxPart;
    int bunchSize = All.BunchSize;
    int maxNodes = MaxNodes;
    integertime ti_Current = All.Ti_Current;
    double errTol2 = All.ErrTolTheta * All.ErrTolTheta;
    
#ifdef DO_NOT_BRACH_IF
    double dxx, dyy, dzz, pdxx, pdyy, pdzz;
#endif
    
#ifdef RT_USE_GRAVTREE
    double mass_stellarlum[N_RT_FREQ_BINS];
    int k_freq; for(k_freq=0;k_freq<N_RT_FREQ_BINS;k_freq++) {mass_stellarlum[k_freq]=0;}
    double dx_stellarlum=0, dy_stellarlum=0, dz_stellarlum=0, sigma_eff=0;
    int valid_gas_particle_for_rt = 0;
#endif
    
    
    
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE)
    double soft=0, pmass;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    double h_p_inv=0, h_p3_inv=0, u_p=0, zeta, ptype_sec=-1, zeta_sec=0;
#endif
#endif
#ifdef EVALPOTENTIAL
    double facpot;
    MyLongDouble pot;
    pot = 0;
#endif
    
    acc_x = 0;
    acc_y = 0;
    acc_z = 0;
    ninteractions = 0;
    nodesinlist = 0;
    
    
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
        pmass = P[target].Mass;
#endif
        aold = All.ErrTolForceAcc * P[target].OldAcc;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL)
        soft = All.ForceSoftening[ptype];
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORGAS)
        if((ptype == 0) && (PPP[target].Hsml>All.ForceSoftening[ptype]))
        {
            soft = PPP[target].Hsml;
            zeta = PPPZ[target].AGS_zeta;
        } else {
            soft = All.ForceSoftening[ptype];
            zeta = 0;
        }
#endif
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
        if(PPP[target].AGS_Hsml > All.ForceSoftening[ptype])
        {
            soft = PPP[target].AGS_Hsml;
            zeta = PPPZ[target].AGS_zeta;
        } else {
            soft = All.ForceSoftening[ptype];
            zeta = 0;
        }
#endif
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
#if defined(RT_USE_GRAVTREE) || defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
        pmass = GravDataGet[target].Mass;
#endif
        ptype = GravDataGet[target].Type;
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE)
        soft = GravDataGet[target].Soft;
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
        zeta = GravDataGet[target].AGS_zeta;
#endif
#endif
    }
    
#if defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS)
    /* quick check if particle has mass: if not, we won't deal with it */
    if(pmass<=0) return 0;
#endif
    
    
    
    
    
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
#ifdef ADAPTIVE_GRAVSOFT_FORALL
    h=soft;
#else
    if(ptype==0) h=soft; else h=All.ForceSoftening[ptype];
#endif
    h_inv = 1.0 / h;
    h3_inv = h_inv * h_inv * h_inv;
#endif
    
    
    
#ifdef RT_USE_GRAVTREE
    if(ptype==0) {if((soft>0)&&(pmass>0)) {valid_gas_particle_for_rt = 1;}}
#endif
    
    
    
    
    if(mode == 0)
    {
        no = maxPart;		/* root node */
    }
    else
    {
        nodesinlist++;
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < maxPart)
            {
                /* the index of the node is the index of the particle */
                if(P[no].Ti_current != ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    drift_particle(no, ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }
                
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
#ifdef PERIODIC
                NEAREST_XYZ(dx,dy,dz,-1);
#endif
                r2 = dx * dx + dy * dy + dz * dz;
                
                mass = P[no].Mass;
                

                
#ifdef RT_USE_GRAVTREE
                if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
                {
                    dx_stellarlum=dx; dy_stellarlum=dy; dz_stellarlum=dz;
                    double lum[N_RT_FREQ_BINS];
                    int active_check = rt_get_source_luminosity(no,sigma_eff,lum);
                    int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {if(active_check) {mass_stellarlum[kf]=lum[kf];} else {mass_stellarlum[kf]=0;}}
                }
#endif
                
                
                
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                /* set secondary softening and zeta term */
                ptype_sec = P[no].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
                if(ptype_sec == 0)
#else
                    if(ptype_sec > -1) /* trigger for all particles */
#endif
                    {
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                        double h_no = PPP[no].AGS_Hsml;
#else
                        double h_no = PPP[no].Hsml;
#endif
                        if(h_no > All.ForceSoftening[P[no].Type])
                        {
                            h_p_inv = 1.0 / h_no;
                            zeta_sec = PPPZ[no].AGS_zeta;
                        }
                        else
                        {
                            h_p_inv = 1.0 / All.ForceSoftening[P[no].Type];
                            zeta_sec = 0;
                        }
                    }
                    else
                    {
                        h_p_inv = 1.0 / All.ForceSoftening[P[no].Type];
                        zeta_sec = 0;
                    }
#else
                h = All.ForceSoftening[ptype];
                if(h < All.ForceSoftening[P[no].Type])
                    h = All.ForceSoftening[P[no].Type];
#endif
                
                
                
                
                
                if(TakeLevel >= 0)
                {
                    LOCK_WORKCOUNT;
#ifdef _OPENMP
#ifdef THREAD_SAFE_COSTS
#pragma omp critical(_workcount_)
#endif
#endif
                    P[no].GravCost[TakeLevel] += 1.0;
                    UNLOCK_WORKCOUNT;
                }
                no = Nextnode[no];
            }
            else			/* we have an  internal node */
            {
                if(no >= maxPart + maxNodes)	/* pseudo particle */
                {
                    if(mode == 0)
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
                
                nop = &Nodes[no];
                
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                
                mass = nop->u.d.mass;
                
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    if(mass)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                
                if(nop->Ti_current != ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    force_drift_node(no, ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }
                
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
                NEAREST_XYZ(dx,dy,dz,-1);
#endif
                r2 = dx * dx + dy * dy + dz * dz;
                

                
#ifdef RT_USE_GRAVTREE
                if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
                {
                    int kf; for(kf=0;kf<N_RT_FREQ_BINS;kf++) {mass_stellarlum[kf] = nop->stellar_lum[kf];}
                    dx_stellarlum = dx; dy_stellarlum = dy; dz_stellarlum = dz;
                }
#endif
                
                
                
                
                
                if(errTol2)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * errTol2)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                else		/* check relative opening criterion */
                {
                    /* force node to open if we are within the gravitational softening length */
#if !(defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE))
                    double soft = All.ForceSoftening[ptype];
#endif
                    if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
#ifdef DO_NOT_BRACH_IF
                    if((mass * nop->len * nop->len > r2 * r2 * aold) |
                       ((pdxx < 0.60 * nop->len) & (pdyy < 0.60 * nop->len) & (pdzz < 0.60 * nop->len)))
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
#else
                    if(mass * nop->len * nop->len > r2 * r2 * aold)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    /* check in addition whether we lie inside the cell */
                    
                    if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
                    {
                        if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
                        {
                            if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
                            {
                                no = nop->u.d.nextnode;
                                continue;
                            }
                        }
                    }
#endif
                }
                
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                /* set secondary softening and zeta term */
                if (nop->maxsoft > 0) h_p_inv = 1.0 / nop->maxsoft; else h_p_inv = 0;
                zeta_sec = 0;
                ptype_sec = -1;
                
                if(h < nop->maxsoft) // compare primary softening to node maximum
                {
                    if(r2 < nop->maxsoft * nop->maxsoft) // inside node maxsoft! continue down tree
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#else
                h = All.ForceSoftening[ptype];
                if(h < nop->maxsoft)
                {
                    h = nop->maxsoft;
                    if(r2 < h * h)
                    {
                        if(maskout_different_softening_flag(nop->u.d.bitflags))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                    }
                }
#endif
                
                if(TakeLevel >= 0)
                {
                    LOCK_WORKCOUNT;
#ifdef _OPENMP
#ifdef THREAD_SAFE_COSTS
#pragma omp critical(_workcount_)
#endif
#endif
                    nop->GravCost += 1.0;
                    UNLOCK_WORKCOUNT;
                }
                
                no = nop->u.d.sibling;	/* ok, node can be used */
                
            }
            
            r = sqrt(r2);
            
            if(r >= h)
            {
                fac = mass / (r2 * r);
#ifdef EVALPOTENTIAL
                facpot = -mass / r;
#endif
            }
            else
            {
#if !defined(ADAPTIVE_GRAVSOFT_FORALL) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
                h_inv = 1.0 / h;
                h3_inv = h_inv * h_inv * h_inv;
#endif
                u = r * h_inv;
                fac = mass * kernel_gravity(u, h_inv, h3_inv, 1);
                
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                // first, appropriately symmetrize the forces between particles //
                if((h_p_inv > 0) && (ptype_sec > -1))
                {
#ifdef HYDRO_SPH
                    if(h_p_inv != h_inv)
                    {
                        h_p3_inv = h_p_inv * h_p_inv * h_p_inv;
                        u_p = r * h_p_inv;
                        fac += mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1);
                        fac *= 0.5;
                    } else {
                        if(zeta_sec != 0)
                        {u_p=u; h_p3_inv=h3_inv;}
                    }
                    // correction only applies to 'shared-kernel' particles: so this needs to check if
                    // these are the same particles for which the kernel lengths are computed
                    if(ags_gravity_kernel_shared_check(ptype, ptype_sec))
                    {
                        double dWdr, wp;
                        if((r>0) && (u<1) && (pmass>0)) // checks that these aren't the same particle
                        {
                            kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
                            fac -= (zeta/pmass) * dWdr / r;   // 0.5 * zeta * omega * dWdr / r;
                        } // if(ptype==0)
                        
                        if(zeta_sec != 0) // secondary is adaptively-softened particle (set above)
                            if(h_p_inv > 0)
                                if((r>0) && (u_p<1) && (pmass>0))
                                {
                                    kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
                                    fac -= (zeta_sec/pmass) * dWdr / r;
                                } // if(zeta_sec != 0)
                    } // if(ptype==ptype_sec)
#else
                    if(h_p_inv < h_inv)
                    {
                        h_p3_inv = h_p_inv * h_p_inv * h_p_inv;
                        u_p = r * h_p_inv;
                        fac = mass * kernel_gravity(u_p, h_p_inv, h_p3_inv, 1);
                    }
                    // correction only applies to 'shared-kernel' particles: so this needs to check if
                    // these are the same particles for which the kernel lengths are computed
                    // (also checks that these aren't the same particle)
#if (defined(MAGNETIC) || defined(COOLING) || defined(GALSF) || defined(FLAG_NOT_IN_PUBLIC_CODE))
                    /* since these modules imply nonstandard cross-particel interactions for certain types, need to limit the correction terms here */
                    if((ptype>0) && (ptype<4) && (ptype_sec>0) && (ptype_sec<4) && (r > 0) && (pmass > 0))
#else
                    if((r > 0) && (pmass > 0))
#endif
                    {
                        if(ags_gravity_kernel_shared_check(ptype, ptype_sec))
                        {
                            double dWdr, wp;
                            if(h_p_inv >= h_inv)
                            {
                                if((zeta != 0) && (u < 1))
                                {
                                    kernel_main(u, h3_inv, h3_inv*h_inv, &wp, &dWdr, 1);
                                    fac -= 2. * (zeta/pmass) * dWdr / sqrt(r2 + 0.0001/(h_inv*h_inv));   // 0.5 * zeta * omega * dWdr / r;
                                }
                            } else {
                                if((zeta_sec != 0) && (u_p < 1)) // secondary is adaptively-softened particle (set above)
                                {
                                    kernel_main(u_p, h_p3_inv, h_p3_inv*h_p_inv, &wp, &dWdr, 1);
                                    fac -= 2. * (zeta_sec/pmass) * dWdr / sqrt(r2 + 0.0001/(h_p_inv*h_p_inv));
                                }
                            }
                        } // if(ptype==ptype_sec)
                    }
#endif
                }
                
#endif // #if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) //
                
#ifdef EVALPOTENTIAL
                facpot = mass * kernel_gravity(u, h_inv, h3_inv, -1);
#endif
            } // closes r < h clause
            
            {
                
                
#ifdef EVALPOTENTIAL
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
                pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
#endif
                
                acc_x += FLT(dx * fac);
                acc_y += FLT(dy * fac);
                acc_z += FLT(dz * fac);
                
#ifdef EVALPOTENTIAL
                pot += FLT(facpot * shortrange_table_potential[tabindex]);
#endif
            } // closes TABINDEX<NTAB
            
            ninteractions++;
            
            
#ifdef RT_USE_GRAVTREE
            if(valid_gas_particle_for_rt)	/* we have a (valid) gas particle as target */
            {
                r2 = dx_stellarlum*dx_stellarlum + dy_stellarlum*dy_stellarlum + dz_stellarlum*dz_stellarlum; r = sqrt(r2);
                if(r >= soft) {fac=1./(r2*r);} else {h_inv=1./soft; h3_inv=h_inv*h_inv*h_inv; u=r*h_inv; fac=kernel_gravity(u,h_inv,h3_inv,1);}
                if((soft>r)&&(soft>0)) fac *= (r2/(soft*soft)); // don't allow cross-section > r2
                
                
            }
#endif // RT_USE_GRAVTREE
            
            
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                {
                    nodesinlist++;
                    no = Nodes[no].u.d.nextnode;	/* open it */
                }
            }
        }
    }
    
    
    
    /* store result at the proper place */
    if(mode == 0)
    {
        P[target].GravAccel[0] = acc_x;
        P[target].GravAccel[1] = acc_y;
        P[target].GravAccel[2] = acc_z;
#ifdef EVALPOTENTIAL
        P[target].Potential = pot;
#endif
    }
    else
    {
        GravDataResult[target].Acc[0] = acc_x;
        GravDataResult[target].Acc[1] = acc_y;
        GravDataResult[target].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
        GravDataResult[target].Potential = pot;
#endif
        *exportflag = nodesinlist;
    }
    
    return ninteractions;
}





#ifdef PERIODIC
/*! This function computes the Ewald correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly
 *  as it was done in Gadget-1.1. Instead, the tree is walked a second
 *  time. This is actually faster because the "Ewald-Treewalk" can use a
 *  different opening criterion than the normal tree walk. In particular,
 *  the Ewald correction is negligible for particles that are very close,
 *  but it is large for particles that are far away (this is quite
 *  different for the normal direct force). So we can here use a different
 *  opening criterion. Sufficient accuracy is usually obtained if the node
 *  length has dropped to a certain fraction ~< 0.25 of the
 *  BoxLength. However, we may only short-cut the interaction list of the
 *  normal full Ewald tree walk if we are sure that the whole node and all
 *  daughter nodes "lie on the same side" of the periodic boundary,
 *  i.e. that the real tree walk would not find a daughter node or particle
 *  that was mapped to a different nearest neighbour position when the tree
 *  walk would be further refined.
 */
int force_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount,
                                        int *exportindex)
{
    struct NODE *nop = 0;
    int no, cost, listindex = 0;
    double dx, dy, dz, mass, r2;
    int signx, signy, signz, nexp;
    int i, j, k, openflag, task;
    double u, v, w;
    double f1, f2, f3, f4, f5, f6, f7, f8;
    MyLongDouble acc_x, acc_y, acc_z;
    double boxsize, boxhalf;
    double pos_x, pos_y, pos_z, aold;
    
    boxsize = All.BoxSize;
    boxhalf = 0.5 * All.BoxSize;
    
    acc_x = 0;
    acc_y = 0;
    acc_z = 0;
    cost = 0;
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        aold = All.ErrTolForceAcc * P[target].OldAcc;
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
    }
    
    if(mode == 0)
    {
        no = All.MaxPart;		/* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)	/* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign */
                if(P[no].Ti_current != All.Ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    drift_particle(no, All.Ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }
                
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
                mass = P[no].Mass;
            }
            else			/* we have an  internal node */
            {
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
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
                                if(Nexport >= All.BunchSize)
                                {
                                    /* out if buffer space. Need to discard work for this particle and interrupt */
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
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        
                        if(exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                nop = &Nodes[no];
                
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    no = nop->u.d.nextnode;
                    continue;
                }
                
                if(nop->Ti_current != All.Ti_Current)
                {
                    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
                    force_drift_node(no, All.Ti_Current);
                    UNLOCK_PARTNODEDRIFT;
                }
                
                mass = nop->u.d.mass;
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
            }
            
            NEAREST_XYZ(dx,dy,dz,-1);
            
            if(no < All.MaxPart)
                no = Nextnode[no];
            else			/* we have an  internal node. Need to check opening criterion */
            {
                openflag = 0;
                r2 = dx * dx + dy * dy + dz * dz;
                if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        openflag = 1;
                    }
                }
                else		/* check relative opening criterion */
                {
                    if(mass * nop->len * nop->len > r2 * r2 * aold)
                    {
                        openflag = 1;
                    }
                    else
                    {
                        if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
                        {
                            if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
                            {
                                if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
                                {
                                    openflag = 1;
                                }
                            }
                        }
                    }
                }
                
                if(openflag)
                {
                    /* now we check if we can avoid opening the cell */
                    
                    u = nop->center[0] - pos_x;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    u = nop->center[1] - pos_y;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    u = nop->center[2] - pos_z;
                    if(u > boxhalf)
                        u -= boxsize;
                    if(u < -boxhalf)
                        u += boxsize;
                    if(fabs(u) > 0.5 * (boxsize - nop->len))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    /* if the cell is too large, we need to refine
                     * it further
                     */
                    if(nop->len > 0.20 * boxsize)
                    {
                        /* cell is too large */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                
                no = nop->u.d.sibling;	/* ok, node can be used */
            }
            
            /* compute the Ewald correction force */
            
            if(dx < 0)
            {
                dx = -dx;
                signx = +1;
            }
            else
                signx = -1;
            if(dy < 0)
            {
                dy = -dy;
                signy = +1;
            }
            else
                signy = -1;
            if(dz < 0)
            {
                dz = -dz;
                signz = +1;
            }
            else
                signz = -1;
            u = dx * fac_intp;
            i = (int) u;
            if(i >= EN)
                i = EN - 1;
            u -= i;
            v = dy * fac_intp;
            j = (int) v;
            if(j >= EN)
                j = EN - 1;
            v -= j;
            w = dz * fac_intp;
            k = (int) w;
            if(k >= EN)
                k = EN - 1;
            w -= k;
            /* compute factors for trilinear interpolation */
            f1 = (1 - u) * (1 - v) * (1 - w);
            f2 = (1 - u) * (1 - v) * (w);
            f3 = (1 - u) * (v) * (1 - w);
            f4 = (1 - u) * (v) * (w);
            f5 = (u) * (1 - v) * (1 - w);
            f6 = (u) * (1 - v) * (w);
            f7 = (u) * (v) * (1 - w);
            f8 = (u) * (v) * (w);
            acc_x += FLT(mass * signx * (fcorrx[i][j][k] * f1 +
                                         fcorrx[i][j][k + 1] * f2 +
                                         fcorrx[i][j + 1][k] * f3 +
                                         fcorrx[i][j + 1][k + 1] * f4 +
                                         fcorrx[i + 1][j][k] * f5 +
                                         fcorrx[i + 1][j][k + 1] * f6 +
                                         fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8));
            acc_y +=
            FLT(mass * signy *
                (fcorry[i][j][k] * f1 + fcorry[i][j][k + 1] * f2 +
                 fcorry[i][j + 1][k] * f3 + fcorry[i][j + 1][k + 1] * f4 + fcorry[i +
                                                                                  1]
                 [j][k] * f5 + fcorry[i + 1][j][k + 1] * f6 + fcorry[i + 1][j +
                                                                            1][k] *
                 f7 + fcorry[i + 1][j + 1][k + 1] * f8));
            acc_z +=
            FLT(mass * signz *
                (fcorrz[i][j][k] * f1 + fcorrz[i][j][k + 1] * f2 +
                 fcorrz[i][j + 1][k] * f3 + fcorrz[i][j + 1][k + 1] * f4 + fcorrz[i +
                                                                                  1]
                 [j][k] * f5 + fcorrz[i + 1][j][k + 1] * f6 + fcorrz[i + 1][j +
                                                                            1][k] *
                 f7 + fcorrz[i + 1][j + 1][k + 1] * f8));
            cost++;
        }
        
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                    no = Nodes[no].u.d.nextnode;	/* open it */
            }
        }
    }
    
    /* add the result at the proper place */
    
    if(mode == 0)
    {
        P[target].GravAccel[0] += acc_x;
        P[target].GravAccel[1] += acc_y;
        P[target].GravAccel[2] += acc_z;
    }
    else
    {
        GravDataResult[target].Acc[0] = acc_x;
        GravDataResult[target].Acc[1] = acc_y;
        GravDataResult[target].Acc[2] = acc_z;
    }
    
    return cost;
}
#endif // #ifdef PERIODIC //




/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
/*! This function also computes the short-range potential when the TreePM
 *  algorithm is used. This potential is the Newtonian potential, modified
 *  by a complementary error function.
 */
int force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
{
    struct NODE *nop = 0;
    MyLongDouble pot;
    int no, ptype, task, nexport_save, listindex = 0;
    double r2, dx, dy, dz, mass, r, u, h, h_inv;
    double pos_x, pos_y, pos_z, aold;
    double fac, dxx, dyy, dzz;
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
    double soft = 0;
#endif
    
    nexport_save = *nexport;
    pot = 0;
    if(mode == 0)
    {
        pos_x = P[target].Pos[0];
        pos_y = P[target].Pos[1];
        pos_z = P[target].Pos[2];
        ptype = P[target].Type;
        aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if((ptype == 0) && (PPP[target].Hsml > All.ForceSoftening[ptype]))
        {
            soft = PPP[target].Hsml;
        } else {
            soft = All.ForceSoftening[ptype];
        }
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        if(PPP[target].AGS_Hsml > All.ForceSoftening[ptype])
        {
            soft = PPP[target].AGS_Hsml;
        } else {
            soft = All.ForceSoftening[ptype];
        }
#endif
    }
    else
    {
        pos_x = GravDataGet[target].Pos[0];
        pos_y = GravDataGet[target].Pos[1];
        pos_z = GravDataGet[target].Pos[2];
        ptype = GravDataGet[target].Type;
        aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
        if(ptype == 0)
            soft = GravDataGet[target].Soft;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORALL
        soft = GravDataGet[target].Soft;
#endif
    }
    
    if(mode == 0)
    {
        no = All.MaxPart;		/* root node */
    }
    else
    {
        no = GravDataGet[target].NodeList[0];
        no = Nodes[no].u.d.nextnode;	/* open it */
    }
    
    while(no >= 0)
    {
        while(no >= 0)
        {
            if(no < All.MaxPart)	/* single particle */
            {
                /* the index of the node is the index of the particle */
                /* observe the sign  */
                if(P[no].Ti_current != All.Ti_Current)
                    drift_particle(no, All.Ti_Current);
                dx = P[no].Pos[0] - pos_x;
                dy = P[no].Pos[1] - pos_y;
                dz = P[no].Pos[2] - pos_z;
                mass = P[no].Mass;
            }
            else
            {
                if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
                {
                    if(mode == 0)
                    {
                        if(Exportflag[task = DomainTask[no - (All.MaxPart + MaxNodes)]] != target)
                        {
                            Exportflag[task] = target;
                            Exportnodecount[task] = NODELISTLENGTH;
                        }
                        
                        if(Exportnodecount[task] == NODELISTLENGTH)
                        {
                            if(*nexport >= All.BunchSize)
                            {
                                *nexport = nexport_save;
                                if(nexport_save == 0)
                                    endrun(13002);	/* in this case, the buffer is too small to process even a single particle */
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
                        DomainNodeIndex[no - (All.MaxPart + MaxNodes)];
                        if(Exportnodecount[task] < NODELISTLENGTH)
                            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                    }
                    no = Nextnode[no - MaxNodes];
                    continue;
                }
                
                nop = &Nodes[no];
                if(mode == 1)
                {
                    if(nop->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
                    {
                        no = -1;
                        continue;
                    }
                }
                
                if(!(nop->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
                {
                    /* open cell */
                    no = nop->u.d.nextnode;
                    continue;
                }
                if(nop->Ti_current != All.Ti_Current)
                    force_drift_node(no, All.Ti_Current);
                mass = nop->u.d.mass;
                dx = nop->u.d.s[0] - pos_x;
                dy = nop->u.d.s[1] - pos_y;
                dz = nop->u.d.s[2] - pos_z;
            }
            
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
            NEAREST_XYZ(dx,dy,dz,-1);
#endif
            r2 = dx * dx + dy * dy + dz * dz;
            if(no < All.MaxPart)
            {
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
                /* set softening */
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                h = soft;
#else
                if(ptype == 0)
                    h = soft;
                else
                    h = All.ForceSoftening[ptype];
#endif
#else
                h = All.ForceSoftening[ptype];
                if(h < All.ForceSoftening[P[no].Type])
                    h = All.ForceSoftening[P[no].Type];
#endif
                no = Nextnode[no];
            }
            else			/* we have an internal node. Need to check opening criterion */
            {
                dxx = nop->center[0] - pos_x;	/* observe the sign ! */
                dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
                dzz = nop->center[2] - pos_z;
                NEAREST_XYZ(dxx,dyy,dzz,-1);
                if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
                {
                    if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
                else		/* check relative opening criterion */
                {
                    
                    /* force node to open if we are within the gravitational softening length */
#if defined(NOGRAVITY) || defined(FLAG_NOT_IN_PUBLIC_CODE) || (!(defined(ADAPTIVE_GRAVSOFT_FORALL) || defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(RT_USE_GRAVTREE)))
                    double soft = All.ForceSoftening[ptype];
#endif
                    if((r2 < (soft+0.6*nop->len)*(soft+0.6*nop->len)) || (r2 < (nop->maxsoft+0.6*nop->len)*(nop->maxsoft+0.6*nop->len)))
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }

#ifdef DO_NOT_BRACH_IF
                    if((mass * nop->len * nop->len > r2 * r2 * aold) |
                       ((fabs(dxx) < 0.60 * nop->len) & (fabs(dyy) < 0.60 * nop->len) & (fabs(dzz) <
                                                                                         0.60 * nop->len)))
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
#else
                    if(mass * nop->len * nop->len > r2 * r2 * aold)
                    {
                        /* open cell */
                        no = nop->u.d.nextnode;
                        continue;
                    }
                    
                    if(fabs(dxx) < 0.60 * nop->len)
                    {
                        if(fabs(dyy) < 0.60 * nop->len)
                        {
                            if(fabs(dzz) < 0.60 * nop->len)
                            {
                                no = nop->u.d.nextnode;
                                continue;
                            }
                        }
                    }
#endif // DO_NOT_BRACH_IF //
                }
                
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL)
#ifdef ADAPTIVE_GRAVSOFT_FORALL
                h = soft;
#else
                if(ptype == 0)
                    h = soft;
                else
                    h = All.ForceSoftening[ptype];
#endif
                
                if(h < nop->maxsoft)
                {
                    //h = nop->maxsoft; // only applies if symmetrizing with MAX(h_i,h_j)
                    if(r2 < nop->maxsoft * nop->maxsoft)
                    {
                        no = nop->u.d.nextnode;
                        continue;
                    }
                }
#else
                h = All.ForceSoftening[ptype];
                if(h < nop->maxsoft)
                {
                    h = nop->maxsoft;
                    if(r2 < h * h)
                    {
                        /* bit-5 signals that there are particles of
                         * different softening in the node
                         */
                        if(maskout_different_softening_flag(nop->u.d.bitflags))
                        {
                            no = nop->u.d.nextnode;
                            continue;
                        }
                    }
                }
#endif // #if defined(ADAPTIVE_GRAVSOFT_FORGAS) || defined(ADAPTIVE_GRAVSOFT_FORALL) //
                no = nop->u.d.sibling;	/* node can be used */
            }
            
            r = sqrt(r2);
            {
                fac = 1;
                if(r >= h)
                    pot += FLT(-fac * mass / r);
                
                else
                {
                    h_inv = 1.0 / h;
                    u = r * h_inv;
                    pot += FLT( fac * mass * kernel_gravity(u, h_inv, 1, -1) );
                }
            }
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
            pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
        }
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                no = GravDataGet[target].NodeList[listindex];
                if(no >= 0)
                    no = Nodes[no].u.d.nextnode;	/* open it */
            }
        }
    }
    
    /* store result at the proper place */
#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
    if(mode == 0)
        P[target].Potential = pot;
    else
        PotDataResult[target].Potential = pot;
#endif
    return 0;
}









/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
    int i;
    size_t bytes;
    double allbytes = 0, allbytes_topleaves = 0;
    double u;
    
    tree_allocated_flag = 1;
    DomainNodeIndex = (int *) mymalloc("DomainNodeIndex", bytes = NTopleaves * sizeof(int));
    allbytes_topleaves += bytes;
    MaxNodes = maxnodes;
    if(!(Nodes_base = (struct NODE *) mymalloc("Nodes_base", bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
        printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
        endrun(3);
    }
    allbytes += bytes;
    if(!
       (Extnodes_base =
        (struct extNODE *) mymalloc("Extnodes_base", bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
        printf("failed to allocate memory for %d tree-extnodes (%g MB).\n",
               MaxNodes, bytes / (1024.0 * 1024.0));
        endrun(3);
    }
    allbytes += bytes;
    Nodes = Nodes_base - All.MaxPart;
    Extnodes = Extnodes_base - All.MaxPart;
    if(!(Nextnode = (int *) mymalloc("Nextnode", bytes = (maxpart + NTopnodes) * sizeof(int))))
    {
        printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n",
               maxpart + NTopnodes, bytes / (1024.0 * 1024.0));
        exit(0);
    }
    allbytes += bytes;
    if(!(Father = (int *) mymalloc("Father", bytes = (maxpart) * sizeof(int))))
    {
        printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
        exit(0);
    }
    allbytes += bytes;
    if(first_flag == 0)
    {
        first_flag = 1;
        if(ThisTask == 0)
            printf
            ("\nAllocated %g MByte for BH-tree, and %g Mbyte for top-leaves.  (presently allocted %g MB)\n\n",
             allbytes / (1024.0 * 1024.0), allbytes_topleaves / (1024.0 * 1024.0),
             AllocatedBytes / (1024.0 * 1024.0));
        for(i = 0; i < NTAB; i++)
        {
            u = 3.0 / NTAB * (i + 0.5);
            shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
            shortrange_table_potential[i] = erfc(u);
        }
    }
}


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
    if(tree_allocated_flag)
    {
        myfree(Father);
        myfree(Nextnode);
        myfree(Extnodes_base);
        myfree(Nodes_base);
        myfree(DomainNodeIndex);
        tree_allocated_flag = 0;
    }
}





/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void dump_particles(void)
{
    FILE *fd;
    char buffer[200];
    int i;
    
    sprintf(buffer, "particles%d.dat", ThisTask);
    fd = fopen(buffer, "w");
    my_fwrite(&NumPart, 1, sizeof(int), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].Pos[0], 3, sizeof(MyFloat), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
    for(i = 0; i < NumPart; i++)
        my_fwrite(&P[i].ID, 1, sizeof(int), fd);
    fclose(fd);
}



#ifdef PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin. These corrections are obtained by Ewald summation. (See
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, the Ewald correction is not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization.  The Ewald summation is done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void ewald_init(void)
{
#ifndef NOGRAVITY
    int i, j, k, beg, len, size, n, task, count;
    double x[3], force[3];
    char buf[200];
    FILE *fd;
    
    if(ThisTask == 0)
    {
        printf("initialize Ewald correction...\n");
        fflush(stdout);
    }
    
#ifdef DOUBLEPRECISION
    sprintf(buf, "ewald_spc_table_%d_dbl.dat", EN);
#else
    sprintf(buf, "ewald_spc_table_%d.dat", EN);
#endif
    if((fd = fopen(buf, "r")))
    {
        if(ThisTask == 0)
        {
            printf("\nreading Ewald tables from file `%s'\n", buf);
            fflush(stdout);
        }
        
        my_fread(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        my_fread(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
        fclose(fd);
    }
    else
    {
        if(ThisTask == 0)
        {
            printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);
            fflush(stdout);
        }
        
        /* ok, let's recompute things. Actually, we do that in parallel. */
        
        size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;
        beg = ThisTask * size;
        len = size;
        if(ThisTask == (NTask - 1))
            len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
        for(i = 0, count = 0; i <= EN; i++)
            for(j = 0; j <= EN; j++)
                for(k = 0; k <= EN; k++)
                {
                    n = (i * (EN + 1) + j) * (EN + 1) + k;
                    if(n >= beg && n < (beg + len))
                    {
                        if(ThisTask == 0)
                        {
                            if((count % (len / 20)) == 0)
                            {
                                printf("%4.1f percent done\n", count / (len / 100.0));
                                fflush(stdout);
                            }
                        }
                        
                        x[0] = 0.5 * ((double) i) / EN;
                        x[1] = 0.5 * ((double) j) / EN;
                        x[2] = 0.5 * ((double) k) / EN;
                        ewald_force(i, j, k, x, force);
                        fcorrx[i][j][k] = force[0];
                        fcorry[i][j][k] = force[1];
                        fcorrz[i][j][k] = force[2];
                        if(i + j + k == 0)
                            potcorr[i][j][k] = 2.8372975;
                        else
                            potcorr[i][j][k] = ewald_psi(x);
                        count++;
                    }
                }
        
        for(task = 0; task < NTask; task++)
        {
            beg = task * size;
            len = size;
            if(task == (NTask - 1))
                len = (EN + 1) * (EN + 1) * (EN + 1) - beg;
            MPI_Bcast(&fcorrx[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&fcorry[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&fcorrz[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
            MPI_Bcast(&potcorr[0][0][beg], len * sizeof(MyFloat), MPI_BYTE, task, MPI_COMM_WORLD);
        }
        
        if(ThisTask == 0)
        {
            printf("\nwriting Ewald tables to file `%s'\n", buf);
            fflush(stdout);
            if((fd = fopen(buf, "w")))
            {
                my_fwrite(&fcorrx[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&fcorry[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&fcorrz[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                my_fwrite(&potcorr[0][0][0], sizeof(MyFloat), (EN + 1) * (EN + 1) * (EN + 1), fd);
                fclose(fd);
            }
        }
    }
    
    fac_intp = 2 * EN / All.BoxSize;
    for(i = 0; i <= EN; i++)
        for(j = 0; j <= EN; j++)
            for(k = 0; k <= EN; k++)
            {
                potcorr[i][j][k] /= All.BoxSize;
                fcorrx[i][j][k] /= All.BoxSize * All.BoxSize;
                fcorry[i][j][k] /= All.BoxSize * All.BoxSize;
                fcorrz[i][j][k] /= All.BoxSize * All.BoxSize;
            }
    
    if(ThisTask == 0)
    {
        printf("initialization of periodic boundaries finished.\n");
        fflush(stdout);
    }
#endif // #ifndef NOGRAVITY
}


/*! This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images. We here use tri-linear
 *  interpolation to get it from the precomputed table, which contains
 *  one octant around the target particle at the origin. The other
 *  octants are obtained from it by exploiting symmetry properties.
 */
double ewald_pot_corr(double dx, double dy, double dz)
{
    int i, j, k;
    double u, v, w;
    double f1, f2, f3, f4, f5, f6, f7, f8;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    if(dz < 0)
        dz = -dz;
    u = dx * fac_intp;
    i = (int) u;
    if(i >= EN)
        i = EN - 1;
    u -= i;
    v = dy * fac_intp;
    j = (int) v;
    if(j >= EN)
        j = EN - 1;
    v -= j;
    w = dz * fac_intp;
    k = (int) w;
    if(k >= EN)
        k = EN - 1;
    w -= k;
    f1 = (1 - u) * (1 - v) * (1 - w);
    f2 = (1 - u) * (1 - v) * (w);
    f3 = (1 - u) * (v) * (1 - w);
    f4 = (1 - u) * (v) * (w);
    f5 = (u) * (1 - v) * (1 - w);
    f6 = (u) * (1 - v) * (w);
    f7 = (u) * (v) * (1 - w);
    f8 = (u) * (v) * (w);
    return potcorr[i][j][k] * f1 +
    potcorr[i][j][k + 1] * f2 +
    potcorr[i][j + 1][k] * f3 +
    potcorr[i][j + 1][k + 1] * f4 +
    potcorr[i + 1][j][k] * f5 +
    potcorr[i + 1][j][k + 1] * f6 + potcorr[i + 1][j + 1][k] * f7 + potcorr[i + 1][j + 1][k + 1] * f8;
}



/*! This function computes the potential correction term by means of Ewald
 *  summation.
 */
double ewald_psi(double x[3])
{
    double alpha, psi;
    double r, sum1, sum2, hdotx;
    double dx[3];
    int i, n[3], h[3], h2;
    
    alpha = 2.0;
    for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
        for(n[1] = -4; n[1] <= 4; n[1]++)
            for(n[2] = -4; n[2] <= 4; n[2]++)
            {
                for(i = 0; i < 3; i++)
                    dx[i] = x[i] - n[i];
                r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
                sum1 += erfc(alpha * r) / r;
            }
    
    for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
        for(h[1] = -4; h[1] <= 4; h[1]++)
            for(h[2] = -4; h[2] <= 4; h[2]++)
            {
                hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
                h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
                if(h2 > 0)
                    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
            }
    
    r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;
    return psi;
}


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
    double alpha, r2;
    double r, val, hdotx, dx[3];
    int i, h[3], n[3], h2;
    
    alpha = 2.0;
    for(i = 0; i < 3; i++)
        force[i] = 0;
    if(iii == 0 && jjj == 0 && kkk == 0)
        return;
    r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    for(i = 0; i < 3; i++)
        force[i] += x[i] / (r2 * sqrt(r2));
    for(n[0] = -4; n[0] <= 4; n[0]++)
        for(n[1] = -4; n[1] <= 4; n[1]++)
            for(n[2] = -4; n[2] <= 4; n[2]++)
            {
                for(i = 0; i < 3; i++)
                    dx[i] = x[i] - n[i];
                r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
                val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);
                for(i = 0; i < 3; i++)
                    force[i] -= dx[i] / (r * r * r) * val;
            }
    
    for(h[0] = -4; h[0] <= 4; h[0]++)
        for(h[1] = -4; h[1] <= 4; h[1]++)
            for(h[2] = -4; h[2] <= 4; h[2]++)
            {
                hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
                h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
                if(h2 > 0)
                {
                    val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);
                    for(i = 0; i < 3; i++)
                        force[i] -= h[i] * val;
                }
            }
}
#endif // #ifdef PERIODIC //


