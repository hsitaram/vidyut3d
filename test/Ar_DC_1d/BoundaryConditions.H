#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace plasmachem_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                      int dir, int sgn,
                      Array4<Real> const& phi,
                      Array4<Real> const& bc_arr,
                      Array4<Real> const& robin_a,
                      Array4<Real> const& robin_b,
                      Array4<Real> const& robin_f,
                      GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                      GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                      GpuArray<Real, AMREX_SPACEDIM> dx,
                      const Real time,
                      ProbParm const& prob_parm,
                               Real Tg, Real Pg)
    {
        IntVect cell_int(i,j,k);
        IntVect ghost_cell(i,j,k);
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        
        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;

        if(sgn == -1) 
        { // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = prob_parm.V;
            bc_arr(ghost_cell)  = prob_parm.V;
        }
        else if(sgn == 1)
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
        }
        else
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
        } 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn, 
                    int spec_id, Array4<Real> const &phi, 
                    Array4<Real> const& bc_arr,
                    Array4<Real> const& robin_a,
                    Array4<Real> const& robin_b,
                    Array4<Real> const& robin_f,
                    GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                    GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                    GpuArray<Real, AMREX_SPACEDIM> dx,
                    const Real time,
                    ProbParm const& prob_parm,
                    Real Tg, Real Pg)

    {
        IntVect cell_int(i,j,k);
        IntVect ghost_cell(i,j,k);
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        
        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;
        
        //default
        robin_a(ghost_cell) = 0.0;
        robin_b(ghost_cell) = 1.0;
        robin_f(ghost_cell) = 0.0;

        if(spec_id==EDN_ID)
        {
            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            amrex::Real D=plasmachem::diffusion_coeff(spec_id, phi(cell_int,ETEMP_ID),
                                          efieldmag,
                                          Tg,Pg);
            
            amrex::Real Te=phi(cell_int,ETEMP_ID);
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            amrex::Real flx = ne*cbar/4.0*Real(sgn);

            bc_arr(ghost_cell)=(-flx/D);
        }

        if(spec_id==ARp_ID)
        {
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
        
            amrex::Real mu = plasmachem::mobility(spec_id,phi(cell_int,ETEMP_ID),
                                           efieldmag,
                                           Tg,Pg);

            amrex::Real D = plasmachem::diffusion_coeff(spec_id,phi(cell_int,ETEMP_ID),
                                           efieldmag,
                                           Tg,Pg);
            
            amrex::Real efield_dotn=phi(cell_int,EFX_ID+dir)*outward_normal[dir];
            
            amrex::Real ion_drift_flux=mu*phi(cell_int,spec_id)*efield_dotn;
            
            //if driftvel is pointing towards the wall
            //then flux=v n - D dndx = v n
            //else  flux = v n - D dndx = 0 
            bc_arr(ghost_cell)=(ion_drift_flux >= 0.0)?(-ion_drift_flux/D):0.0;
            amrex::Print()<<"cell_int, outward_normal, ion flux, bc:"<<cell_int<<"\t"<<
            outward_normal[dir]<<"\t"<<ion_drift_flux<<"\t"<<bc_arr(ghost_cell)<<"\n";
        }

        if(spec_id==ARm_ID)
        {
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            
            amrex::Real D=plasmachem::diffusion_coeff(spec_id, phi(cell_int,ETEMP_ID),
                                          efieldmag,
                                          Tg,Pg);
            
            amrex::Real cbar=sqrt(8.0*K_B*Tg/PI/plasmachem::get_molwt(spec_id));
            amrex::Real flx = phi(i,j,k,spec_id)*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
        }

        if(spec_id==EEN_ID)
        {
            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            amrex::Real Te=phi(cell_int,ETEMP_ID);
            
            amrex::Real D=plasmachem::diffusion_coeff(spec_id, Te,
                                          efieldmag,
                                          Tg,Pg);
            
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            
            amrex::Real flx = (2.0*K_B*Te)*ne*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
        }

    }


    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_vel(int i, int j, int k, int dir,
                     int specid,
                     Array4<Real> const& phi,
                     GpuArray<Real, AMREX_SPACEDIM> dx,
                     const Real time,
                     ProbParm const& prob_parm,
                     amrex::Real Tg,amrex::Real Pg)
    {
        IntVect cellid(i, j, k);
        amrex::Real efieldmag=std::sqrt( std::pow(phi(cellid,EFX_ID),2.0)+
                                         std::pow(phi(cellid,EFY_ID),2.0)+
                                         std::pow(phi(cellid,EFZ_ID),2.0) );
        
        amrex::Real mu = plasmachem::mobility(specid,phi(cellid,ETEMP_ID),
                                              efieldmag,Tg,Pg);

        return(mu*phi(cellid,EFX_ID+dir));
    }
}
#endif
