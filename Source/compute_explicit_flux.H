#ifndef _compute_explicit_flux_H_
#define _compute_explicit_flux_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

#define heps 1.0e-30

//===========================================================================
// WENO 5 JS (Original Jiang GS, Shu CW. Efficient Implementation of Weighted
// ENO Schemes.
//             J. Comput. Phys. 1996; 126(1): 202–228.
//             doi: 10.1006/jcph.1996.0130)
// Stencils i-2, i-1, i, i+1, i+2, i+3
//===========================================================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
weno_reconstruct(amrex::Real uim2,
                 amrex::Real uim1,
                 amrex::Real ui,
                 amrex::Real uip1,
                 amrex::Real uip2,
                 amrex::Real uip3,
                 amrex::Real &umhalf,
                 amrex::Real &uphalf)
{
    constexpr amrex::Real eps = 1e-6;
    amrex::Real vl[3];
    amrex::Real vr[3];
    amrex::Real beta[3];
    amrex::Real alpha[3];
    amrex::Real s[6];

    s[0] = uim2;
    s[1] = uim1;
    s[2] = ui;
    s[3] = uip1;
    s[4] = uip2;
    s[5] = uip3;

    beta[2] = (13.0 / 12.0) * pow(s[0] - 2.0 * s[1] + s[2], 2) +
    0.25 * pow(s[0] - 4.0 * s[1] + 3.0 * s[2], 2);
    beta[1] = (13.0 / 12.0) * pow(s[1] - 2.0 * s[2] + s[3], 2) +
    0.25 * (s[1] - s[3]) * (s[1] - s[3]);
    beta[0] = (13.0 / 12.0) * pow(s[2] - 2.0 * s[3] + s[4], 2) +
    0.25 * pow(3.0 * s[2] - 4.0 * s[3] + s[4], 2);

    beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
    beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
    beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

    alpha[2] = beta[2];
    alpha[1] = 6.0 * beta[1];
    alpha[0] = 3.0 * beta[0];
    amrex::Real alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

    vl[2] = 2.0 * s[0] - 7.0 * s[1] + 11.0 * s[2];
    vl[1] = -s[1] + 5.0 * s[2] + 2.0 * s[3];
    vl[0] = 2.0 * s[2] + 5.0 * s[3] - s[4];

    uphalf = (1.0 / 6.0) * alpha1 *
    (alpha[2] * vl[2] + alpha[1] * vl[1] + alpha[0] * vl[0]);

    beta[2] = (13.0 / 12.0) * pow(s[5] - 2.0 * s[4] + s[3], 2) +
    0.25 * pow(s[5] - 4.0 * s[4] + 3.0 * s[3], 2);
    beta[1] = (13.0 / 12.0) * pow(s[4] - 2.0 * s[3] + s[2], 2) +
    0.25 * (s[4] - s[2]) * (s[4] - s[2]);
    beta[0] = (13.0 / 12.0) * pow(s[3] - 2.0 * s[2] + s[1], 2) +
    0.25 * pow(3.0 * s[3] - 4.0 * s[2] + s[1], 2);

    beta[2] = 1.0 / ((eps + beta[2]) * (eps + beta[2]));
    beta[1] = 1.0 / ((eps + beta[1]) * (eps + beta[1]));
    beta[0] = 1.0 / ((eps + beta[0]) * (eps + beta[0]));

    alpha[2] = beta[2];
    alpha[1] = 6.0 * beta[1];
    alpha[0] = 3.0 * beta[0];
    alpha1 = 1.0 / (alpha[2] + alpha[1] + alpha[0]);

    vr[2] = 11.0 * s[3] - 7.0 * s[4] + 2.0 * s[5];
    vr[1] = -s[4] + 5.0 * s[3] + 2.0 * s[2];
    vr[0] = 2.0 * s[3] + 5.0 * s[2] - s[1];

    umhalf = (1.0 / 6.0) * alpha1 *
    (alpha[2] * vr[2] + alpha[1] * vr[1] + alpha[0] * vr[0]);
}
//====================================================================
//MINBEE vesion of limiter for
//WAF scheme, see page 502, 
//Toro's Riemann solver textbook
AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
amrex::Real get_slope_lim(amrex::Real c,amrex::Real diff_left,
                          amrex::Real diff_mid,
                          amrex::Real diff_right)
{
    amrex::Real r,lim;
    r=0.0;
    r=(c>=0.0)?diff_left/(diff_mid+heps):diff_right/(diff_mid+heps);

    lim=0.0;
    if(r<=0.0)
    {
        lim=1.0;
    }
    else if(r>0.0 && r<=1.0)
    {
        //NOTE: something funky can happen
        //if c or advective CFL is greater than 1
        //But we should be safe with 
        //explicit advection :)
        lim=1.0-(1.0-amrex::Math::abs(c))*r;
    }
    else
    {
        lim=amrex::Math::abs(c);
    }
    return(lim);
}
//See TVD Version of WAF schemes
//Page 501, Toro's Riemann solver texbook
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_secondorder_WAF_flux(amrex::Real vel_L,amrex::Real vel_R,
                                     amrex::Real phi_Lm1,amrex::Real phi_L,
                                     amrex::Real phi_R,amrex::Real phi_Rp1,
                                     amrex::Real dx, amrex::Real dt)
{
    amrex::Real diff_left=phi_L-phi_Lm1;
    amrex::Real diff_mid=phi_R-phi_L;
    amrex::Real diff_right=phi_Rp1-phi_R;
    amrex::Real vel_mid=0.5*(vel_L+vel_R);
    amrex::Real vel_sgn = std::copysign(1.0, vel_mid);
    amrex::Real c=vel_mid*dt/dx;
    amrex::Real fL=phi_L*vel_L;
    amrex::Real fR=phi_R*vel_R;

    amrex::Real lim=get_slope_lim(c,diff_left,diff_mid,diff_right);

    amrex::Real flux=0.5*(fL+fR)-0.5*vel_sgn*lim*(fR-fL);

    return(flux);
}
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real get_firstorder_upwind_flux(amrex::Real vel_L,amrex::Real vel_R,
                                       amrex::Real phi_L,amrex::Real phi_R)
{
    amrex::Real fL=phi_L*vel_L;
    amrex::Real fR=phi_R*vel_R;
    amrex::Real vel_mid=0.5*(vel_L+vel_R);
    amrex::Real vel_sgn = std::copysign(1.0, vel_mid);
    amrex::Real flux = fL *  0.5 * (1.0 + vel_sgn) 
    + fR *  0.5 * (1.0 - vel_sgn);
    return(flux);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void compute_flux(
    int i,
    int j,
    int k,
    int sweepdir,
    int specid,
    Array4<Real> const& phi, // state variable
    const GpuArray<int, AMREX_SPACEDIM> bclo,
    const GpuArray<int, AMREX_SPACEDIM> bchi,
    const GpuArray<int, AMREX_SPACEDIM> domlo,
    const GpuArray<int, AMREX_SPACEDIM> domhi,
    Array4<Real> const& fx,
    amrex::Real Tg,
    amrex::Real Pg,
    amrex::Real time,
    const GpuArray<Real, AMREX_SPACEDIM>& dx,
    amrex::Real dt,
    ProbParm const& prob_parm,int hyp_order)
{

    //|--left1--|--left--|--right--|--right1--|
    //|--(i-2)--|--(i-1)--|--i--|--(i+1)--|
    //remember: face i separates (i-1) and i cells

    IntVect face(i,j,k);
    IntVect left(i, j, k);
    IntVect right(i, j, k);
    left[sweepdir] -= 1;

    Real phi_L = phi(left, specid);
    Real phi_R = phi(right, specid);

    Real vel_L = plasmachem_transport::compute_vel(left[0],left[1],left[2],sweepdir,
                                                   specid,phi,dx,time,prob_parm,Tg,Pg);

    Real vel_R = plasmachem_transport::compute_vel(right[0],right[1],right[2],sweepdir,
                                                   specid,phi,dx,time,prob_parm,Tg,Pg);

    if(hyp_order==1)
    {
        fx(face)=get_firstorder_upwind_flux(vel_L,vel_R,phi_L,phi_R);
    }
    else if(hyp_order==2)
    {
        IntVect left1(i, j, k);
        IntVect right1(i, j, k);
        left1[sweepdir] -= 2;
        right1[sweepdir] += 1;
        amrex::Real phi_Lm1 = phi(left1, specid);
        amrex::Real phi_Rp1 = phi(right1, specid);
        amrex::Real dxmin=std::min(dx[0],std::min(dx[1],dx[2]));

        fx(face)=get_secondorder_WAF_flux(vel_L,vel_R,phi_Lm1,phi_L,
                                          phi_R,phi_Rp1,dxmin,dt);
    }
    else if(hyp_order==5)
    {
        IntVect left1(i, j, k);
        IntVect right1(i, j, k);
        IntVect left2(i, j, k);
        IntVect right2(i, j, k);
        left1[sweepdir] -= 2;
        right1[sweepdir] += 1;
        left2[sweepdir] -= 3;
        right2[sweepdir] += 2;
        
        amrex::Real vel_Lm1 = plasmachem_transport::compute_vel(left1[0],left1[1],left1[2],
                                                                sweepdir,specid,
                                                                phi,dx,time,prob_parm,Tg,Pg);
        
        amrex::Real vel_Lm2 = plasmachem_transport::compute_vel(left2[0],left2[1],left2[2],
                                                                sweepdir,specid,
                                                                phi,dx,time,prob_parm,Tg,Pg);
        
        amrex::Real vel_Rp1 = plasmachem_transport::compute_vel(right1[0],right1[1],right1[2],
                                                                sweepdir,specid,
                                                                phi,dx,time,prob_parm,Tg,Pg);
        
        amrex::Real vel_Rp2 = plasmachem_transport::compute_vel(right2[0],right2[1],right2[2],
                                                                sweepdir,specid,
                                                                phi,dx,time,prob_parm,Tg,Pg);

        amrex::Real ulm2=phi(left2,specid);
        amrex::Real ulm1=phi(left1,specid);
        amrex::Real ul=phi(left,specid);
        amrex::Real ur=phi(right,specid);
        amrex::Real urp1=phi(right1,specid);
        amrex::Real urp2=phi(right2,specid);
        
        amrex::Real flm2=ulm2*vel_Lm2;
        amrex::Real flm1=ulm1*vel_Lm1;
        amrex::Real fl=ul*vel_L;
        amrex::Real fr=ur*vel_R;
        amrex::Real frp1=urp1*vel_Rp1;
        amrex::Real frp2=urp2*vel_Rp2;

        amrex::Real fmhalf=0.0;
        amrex::Real fphalf=0.0;
        amrex::Real umhalf=0.0;
        amrex::Real uphalf=0.0;
        amrex::Real lambda_max=amrex::Math::abs(0.5*(vel_L+vel_R));

        weno_reconstruct(flm2,flm1,fl,fr,frp1,frp2,fmhalf,fphalf);
        weno_reconstruct(ulm2,ulm1,ul,ur,urp1,urp2,umhalf,uphalf);

        fx(face)=0.5*(fphalf+lambda_max*uphalf) +
        0.5*(fmhalf-lambda_max*umhalf);
    }
    else
    {
        amrex::Abort("hyp_order not implemented yet\n");
    }

    if( (face[sweepdir]==domlo[sweepdir]) && (bclo[sweepdir]==IHNEUBC) )
    {
        fx(face)=0.0;
    }
    if( (face[sweepdir]==(domhi[sweepdir]+1)) && (bchi[sweepdir]==IHNEUBC) )
    {
        fx(face)=0.0;
    }
}
#endif
