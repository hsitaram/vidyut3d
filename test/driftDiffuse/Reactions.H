#ifndef _REACTIONS_H_
#define _REACTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>
#include <VarDefines.H>

using namespace amrex;
namespace plasmachem_reactions
{
    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_react_source(int i, int j, int k,
                                  int specid,
                Array4<Real> const& phi,
                Array4<Real> const& reactsource,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time,
                ProbParm const& prob_parm)
        {
            if(specid==S1_ID)
            {
               reactsource(i,j,k) = prob_parm.r1;
            }
            else 
            {
               reactsource(i,j,k) = prob_parm.r2;
            }
        }

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_potential_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& source,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time,
                ProbParm const& prob_parm)
        {
            source(i,j,k)=0.0;
        }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_electron_inelastic_heating(int i, int j, int k,
                     Array4<Real> const& phi,
                     GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                     GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                     GpuArray<Real, AMREX_SPACEDIM> dx,
                     const Real time,
                     ProbParm const& prob_parm,
                     amrex::Real Tg,
                     amrex::Real pres)
    {
        return(0.0);
    }
}
#endif
