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
    void  compute_react_source(int i, int j, int k,
                              Array4<Real> const& phi,
                              Array4<Real> const& rxnsrc,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time,
                              ProbParm const& prob_parm,
                              amrex::Real Tg, amrex::Real Pg)
    {
            //including elec energy
            for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
            {
                rxnsrc(i,j,k,sp)=0.0;
            }
            amrex::Real k_arrh        = 2.5e-12;
            amrex::Real T_arrh       = 278508.0;
            amrex::Real N=Pg/K_B/Tg;
            amrex::Real electemp=phi(i,j,k,ETEMP_ID);

            //amrex::Print()<<"N,electemp,eden="<<N<<"\t"<<
            //electemp<<"\t"<<phi(i,j,k,EDN_ID)<<"\n";
            
            amrex::Real src=k_arrh*std::exp(-T_arrh/electemp)*N*phi(i,j,k,EDN_ID);
            rxnsrc(i,j,k,EDN_ID)=src;
            rxnsrc(i,j,k,I_ID)=src;
            rxnsrc(i,j,k,EEN_ID)=-24.0*ECHARGE*src;
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
}
#endif
