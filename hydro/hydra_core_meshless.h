/* --------------------------------------------------------------------------------- */
/* this is the sub-routine where we actually extrapolate quantities to the faces 
    and set up, then solve the pair-wise Riemann problem for the method */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
/* --------------------------------------------------------------------------------- */
{
//#define DO_HALFSTEP_FOR_MESHLESS_METHODS 1
    
    double s_star_ij,s_i,s_j,v_frame[3],n_unit[3];
    double distance_from_i[3],distance_from_j[3];
    face_area_dot_vel=face_vel_i=face_vel_j=Face_Area_Norm=0;
    
    /* --------------------------------------------------------------------------------- */
    /* define volume elements and interface position */
    /* --------------------------------------------------------------------------------- */
    V_j = P[j].Mass / SphP[j].Density;
    /*
#if defined(HYDRO_MESHLESS_FINITE_VOLUME) || defined(CONSTRAINED_GRADIENT_MHD)
    s_star_ij = 0;
#else
    s_star_ij = 0.5 * kernel.r * (PPP[j].Hsml - local.Hsml) / (local.Hsml + PPP[j].Hsml);
#endif
    */
    s_star_ij = 0;
    /* ------------------------------------------------------------------------------------------------------------------- */
    /* now we're ready to compute the volume integral of the fluxes (or equivalently an 'effective area'/face orientation) */
    /* ------------------------------------------------------------------------------------------------------------------- */
    double wt_i,wt_j;
#ifdef AGGRESSIVE_SLOPE_LIMITERS
    wt_i=V_i; wt_j=V_j;
#else
#ifdef COOLING
    //wt_i=wt_j = 2.*V_i*V_j / (V_i + V_j); // more conservatively, could use DMIN(V_i,V_j), but that is less accurate
    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.25) {wt_i=wt_j=2.*V_i*V_j/(V_i+V_j);} else {wt_i=V_i; wt_j=V_j;}
#else
    //wt_i=wt_j = (V_i*PPP[j].Hsml + V_j*local.Hsml) / (local.Hsml+PPP[j].Hsml); // should these be H, or be -effective sizes- //
    if((fabs(V_i-V_j)/DMIN(V_i,V_j))/NUMDIMS > 1.50) {wt_i=wt_j=(V_i*PPP[j].Hsml+V_j*local.Hsml)/(local.Hsml+PPP[j].Hsml);} else {wt_i=V_i; wt_j=V_j;}
#endif
#endif
    if(SphP[j].ConditionNumber*SphP[j].ConditionNumber > 1.0e12 + cnumcrit2)
    {
        /* the effective gradient matrix is ill-conditioned: for stability, we revert to the "RSPH" EOM */
        Face_Area_Norm = -(wt_i*V_i*kernel.dwk_i + wt_j*V_j*kernel.dwk_j) / kernel.r;
        Face_Area_Norm *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
        Face_Area_Vec[0] = Face_Area_Norm * kernel.dp[0];
        Face_Area_Vec[1] = Face_Area_Norm * kernel.dp[1];
        Face_Area_Vec[2] = Face_Area_Norm * kernel.dp[2];
        Face_Area_Norm = Face_Area_Norm * Face_Area_Norm * r2;
    } else {
        /* the effective gradient matrix is well-conditioned: we can safely use the consistent EOM */
        // note the 'default' formulation from Lanson and Vila takes wt_i=V_i, wt_j=V_j; but this assumes negligible variation in h between particles;
        //      it is more accurate to use a centered wt (centered face area), which we get by linear interpolation //
        for(k=0;k<3;k++)
        {
            Face_Area_Vec[k] = kernel.wk_i * wt_i * (local.NV_T[k][0]*kernel.dp[0] + local.NV_T[k][1]*kernel.dp[1] + local.NV_T[k][2]*kernel.dp[2])
                       + kernel.wk_j * wt_j * (SphP[j].NV_T[k][0]*kernel.dp[0] + SphP[j].NV_T[k][1]*kernel.dp[1] + SphP[j].NV_T[k][2]*kernel.dp[2]);
            Face_Area_Vec[k] *= All.cf_atime*All.cf_atime; /* Face_Area_Norm has units of area, need to convert to physical */
            Face_Area_Norm += Face_Area_Vec[k]*Face_Area_Vec[k];
        }
    }
    if(Face_Area_Norm == 0)
    {
        memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
        Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
    } else {
        
        if((Face_Area_Norm<=0)||(isnan(Face_Area_Norm)))
        {
            printf("PANIC! Face_Area_Norm=%g Mij=%g/%g wk_ij=%g/%g Vij=%g/%g dx/dy/dz=%g/%g/%g NVT=%g/%g/%g NVT_j=%g/%g/%g \n",Face_Area_Norm,local.Mass,P[j].Mass,kernel.wk_i,
                   kernel.wk_j,V_i,V_j,kernel.dp[0],kernel.dp[1],kernel.dp[2],local.NV_T[0][0],local.NV_T[0][1],local.NV_T[0][2],SphP[j].NV_T[0][0],SphP[j].NV_T[0][1],
                   SphP[j].NV_T[0][2]);fflush(stdout);
        }
        Face_Area_Norm = sqrt(Face_Area_Norm);
        for(k=0;k<3;k++) {n_unit[k] = Face_Area_Vec[k] / Face_Area_Norm;}
        

        /* --------------------------------------------------------------------------------- */
        /* extrapolate the conserved quantities to the interaction face between the particles */
        /* first we define some useful variables for the extrapolation */
        /* --------------------------------------------------------------------------------- */
        s_i =  0.5 * kernel.r;
        s_j = -0.5 * kernel.r;
        //(simple up-winding formulation: use if desired instead of time-centered fluxes)//
        //delta_halfstep_i = kernel.sound_i*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_i>s_i) {delta_halfstep_i=s_i;}
        //delta_halfstep_j = kernel.sound_j*0.5*dt_hydrostep*cs_t_to_comoving_x; if(delta_halfstep_j>-s_j) {delta_halfstep_j=-s_j;}
#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing
            certain types of noise and oscillations (but not always!) */
        s_i += 0.5 * dt_hydrostep * (local.Vel[0]*kernel.dp[0] + local.Vel[1]*kernel.dp[1] + local.Vel[2]*kernel.dp[2]) * rinv;
        s_j += 0.5 * dt_hydrostep * (VelPred_j[0]*kernel.dp[0] + VelPred_j[1]*kernel.dp[1] + VelPred_j[2]*kernel.dp[2]) * rinv;
#endif
        s_i = s_star_ij - s_i; //+ delta_halfstep_i; /* projection element for gradients */
        s_j = s_star_ij - s_j; //- delta_halfstep_j;
        distance_from_i[0]=kernel.dp[0]*rinv; distance_from_i[1]=kernel.dp[1]*rinv; distance_from_i[2]=kernel.dp[2]*rinv;
        for(k=0;k<3;k++) {distance_from_j[k] = distance_from_i[k] * s_j; distance_from_i[k] *= s_i;}
        //for(k=0;k<3;k++) {v_frame[k] = 0.5 * (VelPred_j[k] + local.Vel[k]);}
        for(k=0;k<3;k++) {v_frame[k] = rinv * (-s_i*VelPred_j[k] + s_j*local.Vel[k]);} // allows for face to be off-center (to second-order)
        // (note that in the above, the s_i/s_j terms are crossed with the opposing velocity terms: this is because the face is closer to the
        //   particle with the smaller smoothing length; so it's values are slightly up-weighted //
        
        
        /* now we do the reconstruction (second-order reconstruction at the face) */
        reconstruct_face_states(local.Density, local.Gradients.Density, SphP[j].Density, SphP[j].Gradients.Density,
                                distance_from_i, distance_from_j, &Riemann_vec.L.rho, &Riemann_vec.R.rho, 1);
        reconstruct_face_states(local.Pressure, local.Gradients.Pressure, SphP[j].Pressure, SphP[j].Gradients.Pressure,
                                distance_from_i, distance_from_j, &Riemann_vec.L.p, &Riemann_vec.R.p, 1);
#ifdef EOS_GENERAL
        reconstruct_face_states(local.InternalEnergyPred, local.Gradients.InternalEnergy, SphP[j].InternalEnergyPred, SphP[j].Gradients.InternalEnergy,
                                distance_from_i, distance_from_j, &Riemann_vec.L.u, &Riemann_vec.R.u, 1);
        reconstruct_face_states(kernel.sound_i, local.Gradients.SoundSpeed, kernel.sound_j, SphP[j].Gradients.SoundSpeed,
                                distance_from_i, distance_from_j, &Riemann_vec.L.cs, &Riemann_vec.R.cs, 1);
#endif
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.Vel[k]-v_frame[k], local.Gradients.Velocity[k], VelPred_j[k]-v_frame[k], SphP[j].Gradients.Velocity[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.v[k], &Riemann_vec.R.v[k], 1);
        }
#ifdef MAGNETIC
        int slim_mode = 1;
#ifdef CONSTRAINED_GRADIENT_MHD
        if((local.ConditionNumber < 0) || (SphP[j].FlagForConstrainedGradients == 0)) {slim_mode = 1;} else {slim_mode = -1;}
#endif
        for(k=0;k<3;k++)
        {
            reconstruct_face_states(local.BPred[k], local.Gradients.B[k], BPred_j[k], SphP[j].Gradients.B[k],
                                    distance_from_i, distance_from_j, &Riemann_vec.L.B[k], &Riemann_vec.R.B[k], slim_mode);
        }
#ifdef DIVBCLEANING_DEDNER
        reconstruct_face_states(local.PhiPred, local.Gradients.Phi, PhiPred_j, SphP[j].Gradients.Phi,
                                distance_from_i, distance_from_j, &Riemann_vec.L.phi, &Riemann_vec.R.phi, 2);
#endif
#endif

#ifdef DO_HALFSTEP_FOR_MESHLESS_METHODS
        /* advance the faces a half-step forward in time (given our leapfrog scheme, this actually has
            very, very weak effects on the errors. nonetheless it does help a small amount in reducing 
            certain types of noise and oscillations */
        double dt_half = 0.5*dt_hydrostep;
        for(k=0;k<3;k++)
        {
            Riemann_vec.R.rho -= dt_half * local.Density * local.Gradients.Velocity[k][k];
            Riemann_vec.L.rho -= dt_half * SphP[j].Density * SphP[j].Gradients.Velocity[k][k];
            Riemann_vec.R.p -= dt_half * GAMMA * local.Pressure * local.Gradients.Velocity[k][k];
            Riemann_vec.L.p -= dt_half * GAMMA * SphP[j].Pressure * SphP[j].Gradients.Velocity[k][k];
            double dv_l_half = -dt_half * local.Gradients.Pressure[k] / local.Density;
            double dv_r_half = -dt_half * SphP[j].Gradients.Pressure[k] / SphP[j].Density;
            /* // this part doesn't need to be done, because our derivatives are co-moving; but if they are not, use it //
             Riemann_vec.R.rho -= dt_half * (local.Vel[k]-v_frame[k]) * local.Gradients.Density[k];
             Riemann_vec.L.rho -= dt_half * (VelPred_j[k]-v_frame[k]) * SphP[j].Gradients.Density[k];
             Riemann_vec.R.p -= dt_half * (local.Vel[k]-v_frame[k]) * local.Gradients.Pressure[k];
             Riemann_vec.L.p -= dt_half * (VelPred_j[k]-v_frame[k]) * SphP[j].Gradients.Pressure[k];
             for(int kx=0;kx<3;kx++)
             {
                dv_r_half -= dt_half * (local.Vel[kx]-v_frame[kx]) * local.Gradients.Velocity[k][kx];
                dv_l_half -= dt_half * (VelPred_j[kx]-v_frame[kx]) * SphP[j].Gradients.Velocity[k][kx];
             }
            */
            Riemann_vec.R.v[k] += 0.5 * (dv_l_half - dv_r_half);
            Riemann_vec.L.v[k] += 0.5 * (dv_r_half - dv_l_half);
            v_frame[k] += 0.5*(dv_l_half + dv_l_half);
        }
#endif
        
        /* we need the face velocities, dotted into the face vector, for correction back to the lab frame */
        for(k=0;k<3;k++) {face_vel_i+=local.Vel[k]*n_unit[k]; face_vel_j+=VelPred_j[k]*n_unit[k];}
        face_vel_i /= All.cf_atime; face_vel_j /= All.cf_atime;
        face_area_dot_vel = rinv*(-s_i*face_vel_j + s_j*face_vel_i);

        /* also will need approach velocities to determine maximum upwind pressure */
        double v2_approach = 0;
        double vdotr2_phys = kernel.vdotr2;
        if(All.ComovingIntegrationOn) {vdotr2_phys -= All.cf_hubble_a2 * r2;}
        vdotr2_phys *= 1/(kernel.r * All.cf_atime);
        if(vdotr2_phys < 0) {v2_approach = vdotr2_phys*vdotr2_phys;}
        double vdotf2_phys = face_vel_i - face_vel_j; // need to be careful of sign here //
        if(vdotf2_phys < 0) {v2_approach = DMAX( v2_approach , vdotf2_phys*vdotf2_phys );}
        
        double press_i_tot = local.Pressure + local.Density * v2_approach;
        double press_j_tot = SphP[j].Pressure + SphP[j].Density * v2_approach;
#ifdef MAGNETIC
        press_i_tot += 0.5 * kernel.b2_i * fac_magnetic_pressure;
        press_j_tot += 0.5 * kernel.b2_j * fac_magnetic_pressure;
#endif
        double press_tot_limiter;
#ifdef MAGNETIC
        press_tot_limiter = 2.0 * 1.1 * All.cf_a3inv * (press_i_tot + press_j_tot);
#else 
        press_tot_limiter = 1.1 * All.cf_a3inv * DMAX( press_i_tot , press_j_tot );
#endif
#ifdef EOS_GENERAL
        press_tot_limiter *= 2.0;
#endif
#ifdef AGGRESSIVE_SLOPE_LIMITERS
        press_tot_limiter *= 100.0; // large number
#endif
        
        
        /* --------------------------------------------------------------------------------- */
        /* Alright! Now we're actually ready to solve the Riemann problem at the particle interface */
        /* --------------------------------------------------------------------------------- */
        Riemann_solver(Riemann_vec, &Riemann_out, n_unit, press_tot_limiter);
        /* before going on, check to make sure we have a valid Riemann solution */
        if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M))||(Riemann_out.P_M>1.4*press_tot_limiter))
        {
            /* go to a linear reconstruction of P, rho, and v, and re-try */
            Riemann_vec.R.p = local.Pressure; Riemann_vec.L.p = SphP[j].Pressure;
            Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
            for(k=0;k<3;k++) {Riemann_vec.R.v[k]=local.Vel[k]-v_frame[k]; Riemann_vec.L.v[k]=VelPred_j[k]-v_frame[k];}
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
            Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef EOS_GENERAL
            Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
            Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
            Riemann_solver(Riemann_vec, &Riemann_out, n_unit, 1.4*press_tot_limiter);
            if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M))||(Riemann_out.P_M>2.0*press_tot_limiter))
            {
                /* ignore any velocity difference between the particles: this should gaurantee we have a positive pressure! */
                Riemann_vec.R.p = local.Pressure; Riemann_vec.L.p = SphP[j].Pressure;
                Riemann_vec.R.rho = local.Density; Riemann_vec.L.rho = SphP[j].Density;
                for(k=0;k<3;k++) {Riemann_vec.R.v[k]=0; Riemann_vec.L.v[k]=0;}
#ifdef MAGNETIC
                for(k=0;k<3;k++) {Riemann_vec.R.B[k]=local.BPred[k]; Riemann_vec.L.B[k]=BPred_j[k];}
#ifdef DIVBCLEANING_DEDNER
                Riemann_vec.R.phi = local.PhiPred; Riemann_vec.L.phi = PhiPred_j;
#endif
#endif
#ifdef EOS_GENERAL
                Riemann_vec.R.u = local.InternalEnergyPred; Riemann_vec.L.u = SphP[j].InternalEnergyPred;
                Riemann_vec.R.cs = kernel.sound_i; Riemann_vec.L.cs = kernel.sound_j;
#endif
                Riemann_solver(Riemann_vec, &Riemann_out, n_unit, 2.0*press_tot_limiter);
                if((Riemann_out.P_M<0)||(isnan(Riemann_out.P_M)))
                {
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
                    printf("Riemann Solver Failed to Find Positive Pressure!: Pmax=%g PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g H_ij=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g BL=%g/%g/%g BR=%g/%g/%g phiL/R=%g/%g \n",
                           press_tot_limiter,Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,local.Hsml,PPP[j].Hsml,
                           local.Vel[0]-v_frame[0],local.Vel[1]-v_frame[1],local.Vel[2]-v_frame[2],
                           VelPred_j[0]-v_frame[0],VelPred_j[1]-v_frame[1],VelPred_j[2]-v_frame[2],
                           n_unit[0],n_unit[1],n_unit[2],
                           Riemann_vec.L.B[0],Riemann_vec.L.B[1],Riemann_vec.L.B[2],
                           Riemann_vec.R.B[0],Riemann_vec.R.B[1],Riemann_vec.R.B[2],
                           Riemann_vec.L.phi,Riemann_vec.R.phi);
#else
                    printf("Riemann Solver Failed to Find Positive Pressure!: Pmax=%g PL/M/R=%g/%g/%g Mi/j=%g/%g rhoL/R=%g/%g vL=%g/%g/%g vR=%g/%g/%g n_unit=%g/%g/%g \n",
                           press_tot_limiter,Riemann_vec.L.p,Riemann_out.P_M,Riemann_vec.R.p,local.Mass,P[j].Mass,Riemann_vec.L.rho,Riemann_vec.R.rho,
                           Riemann_vec.L.v[0],Riemann_vec.L.v[1],Riemann_vec.L.v[2],
                           Riemann_vec.R.v[0],Riemann_vec.R.v[1],Riemann_vec.R.v[2],n_unit[0],n_unit[1],n_unit[2]);
#endif
                    exit(1234);
                }
            }
        } // closes loop of alternative reconstructions if invalid pressures are found //
        
        
        /* --------------------------------------------------------------------------------- */
        /* Calculate the fluxes (EQUATION OF MOTION) -- all in physical units -- */
        /* --------------------------------------------------------------------------------- */
        if((Riemann_out.P_M>0)&&(!isnan(Riemann_out.P_M)))
        {
            if(All.ComovingIntegrationOn) {for(k=0;k<3;k++) v_frame[k] /= All.cf_atime;}
            
#if defined(HYDRO_MESHLESS_FINITE_MASS) && !defined(MAGNETIC)
            double facenorm_pm = Face_Area_Norm * Riemann_out.P_M;
            for(k=0;k<3;k++) {Fluxes.v[k] = facenorm_pm * n_unit[k];} /* total momentum flux */
            Fluxes.p = facenorm_pm * (Riemann_out.S_M + face_area_dot_vel); // default: total energy flux = v_frame.dot.mom_flux //
            
#ifndef AGGRESSIVE_SLOPE_LIMITERS
            /* for MFM, do the face correction for adiabatic flows here */
            int use_entropic_energy_equation = 0;
            double du_new = 0;
            double SM_over_ceff = fabs(Riemann_out.S_M) / DMIN(kernel.sound_i,kernel.sound_j);
            if(SM_over_ceff < epsilon_entropic_eos_big)
            {
                use_entropic_energy_equation = 1;
                double PdV_fac = Riemann_out.P_M * vdotr2_phys / All.cf_a2inv;
                double PdV_i = kernel.dwk_i * V_i*V_i * local.DhsmlNgbFactor * PdV_fac;
                double PdV_j = kernel.dwk_j * V_j*V_j * PPP[j].DhsmlNgbFactor * PdV_fac;
                du_new = 0.5 * (PdV_i - PdV_j + facenorm_pm * (face_vel_i+face_vel_j));
                // check if, for the (weakly) diffusive case, heat is (correctly) flowing from hot to cold after particle averaging (flux-limit) //
                if(SM_over_ceff > epsilon_entropic_eos_small)
                {
                    double du_old = facenorm_pm * (Riemann_out.S_M + face_area_dot_vel);
                    if(local.Pressure/local.Density > SphP[j].Pressure/SphP[j].Density)
                    {
                        double dtoj = -du_old + facenorm_pm * face_vel_j;
                        if(dtoj > 0) {use_entropic_energy_equation=0;} else {
                            if(dtoj > -du_new+facenorm_pm*face_vel_j) {use_entropic_energy_equation=0;}}
                    } else {
                        double dtoi = du_old - facenorm_pm * face_vel_i;
                        if(dtoi > 0) {use_entropic_energy_equation=0;} else {
                            if(dtoi > du_new-facenorm_pm*face_vel_i) {use_entropic_energy_equation=0;}}
                    }
                }
                // alright, if we've come this far, we need to subtract -off- the thermal energy part of the flux, and replace it //
            }
            if(SphP[j].ConditionNumber*SphP[j].ConditionNumber > cnumcrit2) {use_entropic_energy_equation=1;}
            if(use_entropic_energy_equation) {Fluxes.p = du_new;}
#endif
            
#else
            
            /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
             which we do following Pakmor et al. 2011 */
            for(k=0;k<3;k++)
            {
                Riemann_out.Fluxes.p += v_frame[k] * Riemann_out.Fluxes.v[k];
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
                /* Riemann_out->Fluxes.rho is un-modified */
                Riemann_out.Fluxes.p += (0.5*v_frame[k]*v_frame[k])*Riemann_out.Fluxes.rho;
                Riemann_out.Fluxes.v[k] += v_frame[k] * Riemann_out.Fluxes.rho; /* just boost by frame vel (as we would in non-moving frame) */
#endif
            }
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Riemann_out.Fluxes.B[k] += -v_frame[k] * Riemann_out.B_normal_corrected;} /* v dotted into B along the normal to the face (careful of sign here) */
#endif
            
            /* ok now we can actually apply this to the EOM */
#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
            Fluxes.rho = Face_Area_Norm * Riemann_out.Fluxes.rho;
#endif
            Fluxes.p = Face_Area_Norm * Riemann_out.Fluxes.p; // this is really Dt of --total-- energy, need to subtract KE component for e */
            for(k=0;k<3;k++) {Fluxes.v[k] = Face_Area_Norm * Riemann_out.Fluxes.v[k];} // momentum flux (need to divide by mass) //
#ifdef MAGNETIC
            for(k=0;k<3;k++) {Fluxes.B[k] = Face_Area_Norm * Riemann_out.Fluxes.B[k];} // magnetic flux (B*V) //
            Fluxes.B_normal_corrected = -Riemann_out.B_normal_corrected * Face_Area_Norm;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME)
            //Fluxes.phi = Riemann_out.Fluxes.phi * Face_Area_Norm; // after testing, we now prefer the mass-based fluxes, now //
            if(Fluxes.rho < 0)
            {
                Fluxes.phi = Fluxes.rho * Riemann_vec.R.phi;
            } else {
                Fluxes.phi = Fluxes.rho * Riemann_vec.L.phi; // phi_ij, phi_L/R, or midpoint phi?
            }
#endif
#endif // MAGNETIC

#if defined(HYDRO_MESHLESS_FINITE_MASS) && !defined(AGGRESSIVE_SLOPE_LIMITERS)
            /* for MFM, do the face correction for adiabatic flows here */
            double SM_over_ceff = fabs(Riemann_out.S_M) / DMIN(kernel.sound_i,kernel.sound_j); // for now use sound speed here (more conservative) vs magnetosonic speed //
            /* if SM is sufficiently large, we do nothing to the equations */
            if(SM_over_ceff < epsilon_entropic_eos_big)
            {
                /* ok SM is small, we should use adiabatic equations instead */
#ifdef MAGNETIC
                // convert the face pressure P_M to a gas pressure alone (subtract the magnetic term) //
                for(k=0;k<3;k++) {Riemann_out.P_M -= 0.5*Riemann_out.Face_B[k]*Riemann_out.Face_B[k];}
#endif
                int use_entropic_energy_equation = 1;
                double facenorm_pm = Riemann_out.P_M * Face_Area_Norm;
                double PdV_fac = Riemann_out.P_M * vdotr2_phys / All.cf_a2inv;
                double PdV_i = kernel.dwk_i * V_i*V_i * local.DhsmlNgbFactor * PdV_fac;
                double PdV_j = kernel.dwk_j * V_j*V_j * PPP[j].DhsmlNgbFactor * PdV_fac;
                double du_old = facenorm_pm * (Riemann_out.S_M + face_area_dot_vel);
                double du_new = 0.5 * (PdV_i - PdV_j + facenorm_pm * (face_vel_i+face_vel_j));
                // more detailed check for intermediate cases //
                if(SM_over_ceff > epsilon_entropic_eos_small)
                {
                    if(local.Pressure/local.Density > SphP[j].Pressure/SphP[j].Density)
                    {
                        double dtoj = -du_old + facenorm_pm * face_vel_j;
                        if(dtoj > 0) {use_entropic_energy_equation=0;} else {
                            if(dtoj > -du_new+facenorm_pm*face_vel_j) {use_entropic_energy_equation=0;}}
                    } else {
                        double dtoi = du_old - facenorm_pm * face_vel_i;
                        if(dtoi > 0) {use_entropic_energy_equation=0;} else {
                            if(dtoi > du_new-facenorm_pm*face_vel_i) {use_entropic_energy_equation=0;}}
                    }
                }
                if(SphP[j].ConditionNumber*SphP[j].ConditionNumber > cnumcrit2) {use_entropic_energy_equation=1;}
                // alright, if we've come this far, we need to subtract -off- the thermal energy part of the flux, and replace it //
                if(use_entropic_energy_equation) {Fluxes.p += du_new - du_old;}
            }
#endif // closes MFM check // 

#endif // endif for clause opening full fluxes (mfv or magnetic)
        } else {
            /* nothing but bad riemann solutions found! */
            memset(&Fluxes, 0, sizeof(struct Conserved_var_Riemann));
#ifdef DIVBCLEANING_DEDNER
            Riemann_out.phi_normal_mean=Riemann_out.phi_normal_db=0;
#endif
        }
    } // Face_Area_Norm != 0
}
