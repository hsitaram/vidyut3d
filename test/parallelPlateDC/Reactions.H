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
                              int specid,
                              Array4<Real> const& phi,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time,
                              ProbParm const& prob_parm,
                              amrex::Real Tg, amrex::Real Pg)
    {
        amrex::Real src=0.0;
        if(specid==EDN_ID || specid==I_ID)
        {
            amrex::Real k_arrh        = 2.5e-12;
            amrex::Real T_arrh       = 278508.0;
            amrex::Real N=Pg/K_B/Tg;
            amrex::Real electemp=phi(i,j,k,ETEMP_ID);
            //amrex::Print()<<"N,electemp,eden="<<N<<"\t"<<
            //electemp<<"\t"<<phi(i,j,k,EDN_ID)<<"\n";
            src=k_arrh*std::exp(-T_arrh/electemp)*N*phi(i,j,k,EDN_ID);
        }
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
        source(i,j,k)=(phi(i,j,k,EDN_ID)-phi(i,j,k,I_ID))
        *(ECHARGE/EPS0)*prob_parm.ns/prob_parm.phis;
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
                                                   amrex::Real Pg)
    {
        amrex::Real ri=0.0;
        amrex::Real Eiri=0.0;
        
        amrex::Real Ei           = -15.578*ECHARGE; //converting eV to J
        amrex::Real k_arrh       = 2.5e-12;
        amrex::Real T_arrh       = 278508.0;
        amrex::Real N=Pg/K_B/Tg;
        amrex::Real electemp=phi(i,j,k,ETEMP_ID);
        ri=k_arrh*std::exp(-T_arrh/electemp)*N*phi(i,j,k,EDN_ID);
        Eiri=Ei*ri;
        return(Eiri);
    }
}
#endif
