#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <Chemistry.H>
#include <ProbParm.H>

using namespace amrex;

AMREX_INLINE
void amrex_probinit (ProbParm & h_prob_parm, ProbParm & d_prob_parm)
{
    // read problem specific parmparse parameters here
    amrex::ParmParse pp("prob");
    pp.query("loc_x", h_prob_parm.loc_x );
    pp.query("loc_y", h_prob_parm.loc_y );
    pp.query("loc_z", h_prob_parm.loc_z );
    pp.query("rad", h_prob_parm.rad );
    pp.query("ampl", h_prob_parm.ampl );
    pp.query("bg", h_prob_parm.bg );
    pp.query("V1", h_prob_parm.V1 );
    pp.query("V2", h_prob_parm.V2 );
    pp.query("dcoeff", h_prob_parm.dcoeff );

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#else
    std::memcpy(&d_prob_parm, &h_prob_parm, sizeof(ProbParm));
#endif

}

AMREX_GPU_DEVICE
AMREX_INLINE
void initdomaindata(Box const& bx, Array4<Real> const& phi, 
        GeometryData const& geomdata,
        ProbParm *prob_parm)
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);    

    const Real* AMREX_RESTRICT prob_lo = geomdata.ProbLo();
    const Real* AMREX_RESTRICT prob_hi = geomdata.ProbHi();
    const Real* AMREX_RESTRICT dx      = geomdata.CellSize();
    
    for (int k = lo.z; k <= hi.z; ++k) 
    {
        for (int j = lo.y; j <= hi.y; ++j) 
        {
            
            for (int i = lo.x; i <= hi.x; ++i) 
            {
                Real z = prob_lo[2] + (0.5+k) * dx[2];
                Real y = prob_lo[1] + (0.5+j) * dx[1];
                Real x = prob_lo[0] + (0.5+i) * dx[0];

                Real dist2 = std::pow(x-prob_parm->loc_x,2.0)
                      + std::pow(y-prob_parm->loc_y,2.0)
                      + std::pow(z-prob_parm->loc_z,2.0);
                
                phi(i,j,k,I_ID)  = prob_parm->bg + prob_parm->ampl*exp(-dist2/prob_parm->rad/prob_parm->rad);
                phi(i,j,k,E_ID)  = prob_parm->bg;
                phi(i,j,k,EFX_ID) = 0.0;
                phi(i,j,k,EFY_ID) = 0.0;
                phi(i,j,k,EFZ_ID) = 0.0;
                phi(i,j,k,POT_ID) = 0.0;
            }
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void externalbc(const amrex::Real x[AMREX_SPACEDIM],
        const amrex::Real s_int[],
        amrex::Real s_ext[],
        const int idir,
        const int sgn,
        const amrex::Real time,
        amrex::GeometryData const& geomdata)
{
    //default to extrapolation
    for(int c=0;c<NVAR;c++)
    {
        s_ext[c]=s_int[c];
    }
}
#endif
