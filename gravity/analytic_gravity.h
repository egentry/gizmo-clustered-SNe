/*! \file analytic_gravity.h
 *  \brief externally-specified (analytic) gravity goes here
 *
 *  This file contains supplemental code if you want to add an 
 *   -analytic- potential or gravitational force in the code, 
 *   rather than solely relying on the calculated self-gravity
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

void add_analytic_gravitational_forces(void);
void GravAccel_StaticPlummerSphere(void);
void GravAccel_StaticHernquist(void);
void GravAccel_StaticIsothermalSphere(void);
void GravAccel_KeplerianOrbit(void);
void GravAccel_KeplerianTestProblem(void);
void GravAccel_GrowingDiskPotential(void);
void GravAccel_StaticNFW(void);
void GravAccel_RayleighTaylorTest(void);
void GravAccel_ShearingSheet(void);
void GravAccel_PaczynskyWiita(void);


/* master routine which decides which (if any) analytic gravitational forces are applied */
void add_analytic_gravitational_forces()
{
#ifdef ANALYTIC_GRAVITY
#ifdef SHEARING_BOX
    GravAccel_ShearingSheet();            // adds coriolis and centrifugal terms for shearing-sheet approximation
#endif
    //GravAccel_RayleighTaylorTest();     // vertical potential for RT tests
    //GravAccel_StaticPlummerSphere();    // plummer sphere
    //GravAccel_StaticHernquist();        // hernquist sphere
    //GravAccel_StaticIsothermalSphere(); // singular or cored isothermal sphere
    //GravAccel_KeplerianOrbit();         // keplerian disk
    //GravAccel_KeplerianTestProblem();   // keplerian disk with boundaries for test problem
    //GravAccel_GrowingDiskPotential();   // time-dependent (adiabatically growing) disk
    //GravAccel_StaticNFW();              // spherical NFW profile
    //GravAccel_PaczynskyWiita();         // Paczynsky-Wiita pseudo-Newtonian potential
#endif
}


/* adds coriolis and centrifugal terms for shearing-sheet approximation */
void GravAccel_ShearingSheet()
{
#ifdef SHEARING_BOX
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* centrifugal force term (depends on distance from box center) */
        P[i].GravAccel[0] += 2.*(P[i].Pos[0]-boxHalf_X) * SHEARING_BOX_Q*SHEARING_BOX_OMEGA_BOX_CENTER*SHEARING_BOX_OMEGA_BOX_CENTER;
        /* coriolis force terms */
        double vp=0;
        if(P[i].Type==0) {vp=SphP[i].VelPred[SHEARING_BOX_PHI_COORDINATE];} else {vp=P[i].Vel[SHEARING_BOX_PHI_COORDINATE];}
        P[i].GravAccel[0] += 2.*vp * SHEARING_BOX_OMEGA_BOX_CENTER;
        if(P[i].Type==0) {vp=SphP[i].VelPred[0];} else {vp=P[i].Vel[0];}
        P[i].GravAccel[SHEARING_BOX_PHI_COORDINATE] -= 2.*vp * SHEARING_BOX_OMEGA_BOX_CENTER;
#if (SHEARING_BOX==4)
        /* add vertical gravity to the force law */
        P[i].GravAccel[2] -= SHEARING_BOX_OMEGA_BOX_CENTER * SHEARING_BOX_OMEGA_BOX_CENTER * (P[i].Pos[2]-boxHalf_Z);
#endif
    }
#endif
}



/* constant vertical acceleration for Rayleigh-Taylor test problem */
void GravAccel_RayleighTaylorTest()
{
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        /* zero out the gravity first (since this test doesn't use self-gravity) */
        P[i].GravAccel[0]=P[i].GravAccel[1]=P[i].GravAccel[2]=0;
        /* now add the constant vertical field */
        if(P[i].ID != 0) {P[i].GravAccel[1]=-0.5;}
    }
}



/* static unit Plummer Sphere (G=M=a=1) */
void GravAccel_StaticPlummerSphere()
{
    int i,l; double r;
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
        for(l = 0; l < 3; l++)
            P[i].GravAccel[l] += -P[i].Pos[l] / pow(r * r + 1, 1.5);
        
    }
}



/* static Hernquist Profile (parameters specified in the routine below) */
void GravAccel_StaticHernquist()
{
    double HQ_M200 = 95.2401;
    double HQ_C = 9.0;
    double HQ_DARKFRACTION = 0.9;

    double r, m, a; int i,l;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        r = sqrt(P[i].Pos[0] * P[i].Pos[0] + P[i].Pos[1] * P[i].Pos[1] + P[i].Pos[2] * P[i].Pos[2]);
        
        a = pow(All.G * HQ_M200 / (100 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits), 1.0 / 3) / HQ_C *
        sqrt(2 * (log(1 + HQ_C) - HQ_C / (1 + HQ_C)));
        
        m = HQ_M200 * pow(r / (r + a), 2) * HQ_DARKFRACTION;
        if(r > 0)
        {
            for(l = 0; l < 3; l++)
                P[i].GravAccel[l] += -All.G * m * P[i].Pos[l] / (r * r * r);
            
        }
    }
}



/* static Isothermal Sphere Profile (parameters specified in the routine below) */
void GravAccel_StaticIsothermalSphere()
{
    double ISO_M200=95.21;
    double ISO_R200=160.0;
    double ISO_Eps=0.1;
    double ISO_FRACTION=0.9;
    double r, m, dx, dy, dz; int i;
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dx = P[i].Pos[0]; dy = P[i].Pos[1]; dz = P[i].Pos[2];
        r = sqrt(dx * dx + dy * dy + dz * dz);
        
        if(r > ISO_R200)
            m = ISO_M200;
        else
            m = ISO_M200 * r / ISO_R200;
        
        m *= ISO_FRACTION;
        if(r > 0)
        {
            P[i].GravAccel[0] += -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
            P[i].GravAccel[1] += -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
            P[i].GravAccel[2] += -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
        }
    }
}


/* time-dependent potential of an adiabatically-growing disk */
void GravAccel_GrowingDiskPotential()
{
#ifdef GROWING_DISK_POTENTIAL
    double mdisk, dx, dy, dz, r, z, aR, az; int i;
    /* currently ifdef'd out because these routines need to be supplied externally */
    growing_disk_init();
    mdisk = get_disk_mass(All.Time);
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dx = P[i].Pos[0]; dy = P[i].Pos[1]; dz = P[i].Pos[2];
        r = sqrt(dx * dx + dy * dy); z = fabs(dz);
        
        get_disk_forces(r, z, &aR, &az);
        
        aR *= mdisk;
        az *= mdisk;
        
        if(r > 0)
        {
            P[i].GravAccel[0] += -dx / r * aR;
            P[i].GravAccel[1] += -dy / r * aR;
            P[i].GravAccel[2] += -dz / z * az;
        }
    }
#endif
}


/* Keplerian forces (G=M=1): useful for orbit, MRI, planetary disk problems */
void GravAccel_KeplerianOrbit()
{
    double x00,y00;
    x00=y00=0;
#if defined(PERIODIC)
    x00=boxHalf_X; y00=boxHalf_Y;
#endif
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].GravAccel[0] = -(P[i].Pos[0]-x00)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
        P[i].GravAccel[1] = -(P[i].Pos[1]-y00)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
        P[i].GravAccel[2] = 0;
    }
}





/* Keplerian forces (G=M=1): this is a specific (bounded and softened) version 
 used just for the Keplerian disk test problem */
void GravAccel_KeplerianTestProblem()
{
    double x00=0;//boxHalf_X;
    double y00=0;//boxHalf_Y;
    x00=4.0;
    y00=4.0;
    
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double r = pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),0.5);
        if((r > 0.35)&(r < 2.1))
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r <= 0.35)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*pow(r/0.35,2)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*pow(r/0.35,2)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            
            P[i].GravAccel[0] += +(P[i].Pos[0]-x00)*(0.35-r)/0.35
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] += +(P[i].Pos[1]-y00)*(0.35-r)/0.35
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
        if(r >= 2.1)
        {
            P[i].GravAccel[0] = -(P[i].Pos[0]-x00)*(1+(r-2.1)/0.1)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[1] = -(P[i].Pos[1]-y00)*(1+(r-2.1)/0.1)
            / pow(pow(P[i].Pos[1]-y00,2.)+pow(P[i].Pos[0]-x00,2.),1.5) ;
            P[i].GravAccel[2] = 0;
        }
    }
}





/* static NFW potential */
void GravAccel_StaticNFW()
{
    double NFW_C=12;
    double NFW_M200=100.0;
    double NFW_Eps=0.01;
    double NFW_DARKFRACTION=0.87;
    double NFW_BOXCENTERED;
    NFW_BOXCENTERED=1;

    /* convert units */
    double R200 = pow(NFW_M200 * All.G / (100 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits), 1.0 / 3);
    double Rs = R200 / NFW_C;
    double Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
    double RhoCrit = 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    double V200 = 10 * All.Hubble_H0_CodeUnits * R200;
    
    double r0, R, r, m, dx, dy, dz, fac; int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        dx = P[i].Pos[0];
        dy = P[i].Pos[1];
        dz = P[i].Pos[2];
#if defined(PERIODIC)
        if(NFW_BOXCENTERED)
        {
            dx = P[i].Pos[0] - boxHalf_X;
            dy = P[i].Pos[1] - boxHalf_Y;
            dz = P[i].Pos[2] - boxHalf_Z;
        }
#endif
        r0 = sqrt(dx * dx + dy * dy + dz * dz);

        /* function to get enclosed mass(<r) for NFW: */
        /* Eps is in units of Rs !!!! :: use unsoftened NFW if NFW_Eps=0 */
        R = r0;
        if(NFW_Eps > 0.0)
            if(R > Rs * NFW_C)
                R = Rs * NFW_C;
        
        fac=1.0;
        if(NFW_Eps > 0.0)
        {
            m = fac * 4 * M_PI * RhoCrit * Dc * (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
              NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
              NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs)))
              / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
        }
        else /* analytic NFW */
        {
            m = fac * 4 * M_PI * RhoCrit * Dc *
            (-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
        }
        fac = V200 * V200 * V200 / (10 * All.G * All.Hubble_H0_CodeUnits) / m;
        if(NFW_Eps > 0.0)
        {
            m = fac * 4 * M_PI * RhoCrit * Dc * (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) +
              NFW_Eps * NFW_Eps * log(NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1)) + (Rs * Rs * Rs * (Rs -
              NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) + NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs)))
              / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
        }
        else /* analytic NFW */
        {
            m = fac * 4 * M_PI * RhoCrit * Dc *
            (-(Rs * Rs * Rs * (1 + log(Rs))) + Rs * Rs * Rs * (Rs + (R + Rs) * log(R + Rs)) / (R + Rs));
        }
        m *= NFW_DARKFRACTION; r=r0;

        if(r > 0)
        {
            P[i].GravAccel[0] += -All.G * m * dx / (r * r * r);
            P[i].GravAccel[1] += -All.G * m * dy / (r * r * r);
            P[i].GravAccel[2] += -All.G * m * dz / (r * r * r);
            
        } // if(r > 0) //
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //
}




/* Paczysnky Wiita pseudo-Newtonian potential, G = M_sol = c = 1 */
void GravAccel_PaczynskyWiita()
{
    double PACZYNSKY_WIITA_MASS = 1.0; // Mass to use for the Paczynksy-Wiita analytic gravity pseudo-Newtonian potential (in solar masses)
    double r_g = 2*PACZYNSKY_WIITA_MASS;
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        double r = sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
        if(r > r_g)
        {
            double q = PACZYNSKY_WIITA_MASS/((r - r_g)*(r - r_g));
            P[i].GravAccel[0] = - q * P[i].Pos[0]/r;
            P[i].GravAccel[1] = - q * P[i].Pos[1]/r;
            P[i].GravAccel[2] = - q * P[i].Pos[2]/r;
        }
    }
}

#ifdef PARTICLE_EXCISION
void apply_excision(void)
{
    double EXCISION_MASS = 0; // mass of the excised object. Used to move the excision boundary so as to capture bound objects. If zero the excision boundary will not move
    double EXCISION_INIT_RADIUS = 0; // initial excision radius
    double EXCISION_ETA = 1; // remove particles with radius < EXCISION_ETA R_excision
    double excision_radius = EXCISION_ETA * pow(EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS*EXCISION_INIT_RADIUS +
                                                3.*sqrt(2. * All.G * EXCISION_MASS) * pow(EXCISION_INIT_RADIUS, 3./2.) * All.Time +
                                                9./2. * All.G * EXCISION_MASS * All.Time*All.Time, 1./3.);
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0)
        {
            double r = sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
            if(r < excision_radius) P[i].Mass = 0;
        }
    }
}
#endif

