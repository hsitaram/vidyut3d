#ifndef _compute_flux_3d_H_
#define _compute_flux_3d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <Constants.H>
#include <Reactions.H>

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
        Array4<Real> const& vel,
        Array4<Real> const& fx,
        const GpuArray<Real, AMREX_SPACEDIM>& dx,
        ProbParm const& prob_parm)
{

    //|--left1--|--left--|--right--|--right1--|
    //|--(i-2)--|--(i-1)--|--i--|--(i+1)--|
    //remember: face i separates (i-1) and i cells

    IntVect left(i, j, k);
    IntVect left1(i, j, k);
    IntVect right(i, j, k);
    IntVect right1(i, j, k);

    left[sweepdir] -= 1;
    left1[sweepdir] -= 2;
    right1[sweepdir] += 1;

    // find face vel
    Real smallval = 1e-10;
    Real dtr = phi(right, specid) - phi(left, specid);
    if (!(std::abs(dtr) > 0.0))
    {
        dtr = dtr + smallval;
    }
    Real r_left = (phi(left, specid) - phi(left1, specid)) / dtr;
    Real lim_left = std::max(0.0, std::min(1.0, r_left));

    dtr = phi(right1, specid) - phi(right, specid);
    if (!(std::abs(dtr) > 0.0))
    {
        dtr = dtr + smallval;
    }
    Real r_right = (phi(right, specid) - phi(left, specid)) / dtr;
    Real lim_right = std::max(0.0, std::min(1.0, r_right));

    //Real phi_L = phi(left, specid) + 0.5 * lim_left * (phi(left, specid) - phi(left1, specid));
    //Real phi_R = phi(right, specid) - 0.5 * lim_right * (phi(right1, specid) - phi(right, specid));
    
    Real phi_L = phi(left, specid);
    Real phi_R = phi(right, specid);

    Real vel_mid = vel(i, j, k);
    fx(i, j, k) = phi_L * 0.5 * (vel_mid + fabs(vel_mid)) + phi_R * 0.5 * (vel_mid - fabs(vel_mid));
}
#endif
