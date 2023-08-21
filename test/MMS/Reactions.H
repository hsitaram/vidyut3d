#ifndef _REACTIONS_H_
#define _REACTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

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
        int bigdir=0;
        amrex::Real xlen=prob_hi[0]-prob_lo[0];
        amrex::Real ylen=prob_hi[1]-prob_lo[1];
        amrex::Real zlen=prob_hi[2]-prob_lo[2];
        
        amrex::Real maxlen=std::max(xlen,std::max(ylen,zlen));
       
        if(maxlen==xlen) bigdir=0;
        if(maxlen==ylen) bigdir=1;
        if(maxlen==zlen) bigdir=2;
       
        IntVect iv(i, j, k);
        amrex::Real x=prob_lo[bigdir]+(iv[bigdir]+0.5)*dx[bigdir];
        if(specid==S1_ID)
        {
                reactsource(i,j,k) = 5.0*std::pow(x,4.0)/3.0-x/6.0-2.0;
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

        source(i,j,k)=phi(i,j,k,S1_ID);
    }
}
#endif
