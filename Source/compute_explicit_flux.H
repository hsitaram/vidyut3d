#ifndef _compute_explicit_flux_H_
#define _compute_explicit_flux_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

#define heps 1.0e-30

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
    if(c>=0.0)
    {
        r=diff_left/(diff_mid+heps);
    }
    else
    {
        r=diff_right/(diff_mid+heps);
    }

    lim=0.0;
    if(r<=0.0)
    {
        lim=1.0;
    }
    else if(r>0.0 && r<=1.0)
    {
        lim=1.0-(1.0-amrex::Math::abs(c))*r;
    }
    else
    {
        lim=amrex::Math::abs(c);
    }

    return(lim);
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
    IntVect left1(i, j, k);
    IntVect right(i, j, k);
    IntVect right1(i, j, k);

    left[sweepdir] -= 1;

    // find face vel
    Real smallval = 1e-10;

    Real phi_L = phi(left, specid);
    Real phi_R = phi(right, specid);

    Real vel_L = plasmachem_transport::compute_vel(left[0],left[1],left[2],sweepdir,
                                                   specid,phi,dx,time,prob_parm,Tg,Pg);

    Real vel_R = plasmachem_transport::compute_vel(right[0],right[1],right[2],sweepdir,
                                                   specid,phi,dx,time,prob_parm,Tg,Pg);

    Real vel_mid=0.5*(vel_L+vel_R);

    Real fL=phi_L*vel_L;
    Real fR=phi_R*vel_R;
    Real vel_sgn = std::copysign(1.0, vel_mid);

    if(hyp_order==1)
    {

        fx(face) = fL *  0.5 * (1.0 + vel_sgn) 
        + fR *  0.5 * (1.0 - vel_sgn);
    }
    else if(hyp_order==2)
    {
        left1[sweepdir] -= 2;
        right1[sweepdir] += 1;
        Real phi_Lm1 = phi(left1, specid);
        Real phi_Rp1 = phi(right1, specid);
        Real dxmin=std::min(dx[0],std::min(dx[1],dx[2]));

        Real diff_left=phi_L-phi_Lm1;
        Real diff_mid=phi_R-phi_L;
        Real diff_right=phi_Rp1-phi_R;
        Real c=vel_mid*dt/dxmin;

        Real lim=get_slope_lim(c,diff_left,diff_mid,diff_right);

        //See TVD Version of WAF schemes
        //Page 501, Toro's Riemann solver texbook
        fx(face)=0.5*(fL+fR)-0.5*vel_sgn*lim*(fR-fL);
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
