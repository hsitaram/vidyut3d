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
    amrex::Real compute_react_source(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& rxnsrc,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time,
                              ProbParm const& prob_parm,
                              amrex::Real Tg,amrex::Real Pg)
    {

        //Bagheri, Behnaz, et al. "Comparison of six simulation codes for 
        //positive streamers in air." 
        //Plasma Sources Science and Technology 27.9 (2018): 095002.        
        
        //including elec energy
        for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
        {
            rxnsrc(i,j,k,sp)=0.0;
        }

        amrex::Real src=0.0;
        Real elecfield_mag=std::sqrt(std::pow(phi(i,j,k,EFX_ID),2.0) +
                                     std::pow(phi(i,j,k,EFY_ID),2.0) +
                                     std::pow(phi(i,j,k,EFZ_ID),2.0) )*prob_parm.phis;

        Real alpha = (1.1944e6 + 4.3666e26/std::pow(elecfield_mag,3.0))
        *std::exp(-2.73e7/elecfield_mag);
        Real eta = 340.75;
        Real alpha_bar = alpha-eta;
        Real mu_e = 2.3987 * std::pow(elecfield_mag,-0.26);

        Real react_term=alpha_bar*mu_e*elecfield_mag*phi(i,j,k,EDN_ID);

        src=react_term/prob_parm.ns;
        rxnsrc(i,j,k,EDN_ID)=src;
        rxnsrc(i,j,k,I_ID)=src;
        return(src);
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
        source(i,j,k)=(phi(i,j,k,EDN_ID)-phi(i,j,k,I_ID))*(ECHARGE/EPS0)
        *prob_parm.ns/prob_parm.phis;
    }
}
#endif
