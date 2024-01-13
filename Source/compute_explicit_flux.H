#ifndef _compute_explicit_flux_H_
#define _compute_explicit_flux_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <UnivConstants.H>
#include <VarDefines.H>

using namespace amrex;

#define heps 1.0e-30

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
        ProbParm const& prob_parm)
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
    left1[sweepdir] -= 2;
    right1[sweepdir] += 1;

    // find face vel
    Real smallval = 1e-10;
    
    Real phi_L = phi(left, specid);
    Real phi_R = phi(right, specid);

    Real vel_L = plasmachem_transport::compute_vel(left[0],left[1],left[2],sweepdir,
                     specid,phi,dx,time,prob_parm,Tg,Pg);

    Real vel_R = plasmachem_transport::compute_vel(right[0],right[1],right[2],sweepdir,
                     specid,phi,dx,time,prob_parm,Tg,Pg);

    Real vel_mid=0.5*(vel_L+vel_R);
    Real vel_sgn = std::copysign(1.0, vel_mid);

    fx(face) = phi_L * vel_L *  0.5 * (vel_sgn + fabs(vel_sgn)) 
                + phi_R * vel_R *  0.5 * (vel_sgn - fabs(vel_sgn));
    
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
