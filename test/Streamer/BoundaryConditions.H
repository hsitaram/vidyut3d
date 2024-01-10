#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace plasmachem_transport
{
    // AMREX_GPU_DEVICE AMREX_INLINE
    // void potential_bc(int i, int j, int k,
    //                   int dir, int sgn,
    //                   Array4<Real> const& phi,
    //                   Array4<Real> const& bc_arr,
    //                   Array4<Real> const& robin_a,
    //                   Array4<Real> const& robin_b,
    //                   Array4<Real> const& robin_f,
    //                   GpuArray<Real, AMREX_SPACEDIM> prob_lo,
    //                   GpuArray<Real, AMREX_SPACEDIM> prob_hi,
    //                   GpuArray<Real, AMREX_SPACEDIM> dx,
    //                   const Real time,
    //                   ProbParm const& prob_parm,
    //                   amrex::Real Tg,amrex::Real Pg)
    // {
    //     IntVect cell_int(i,j,k);
    //     IntVect ghost_cell(i,j,k);
    //     amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
    //     outward_normal[dir]=sgn;
    //     
    //     int gcell_adjust,cell_adjust;

    //     //ghost_cell is one behind
    //     gcell_adjust= (sgn==-1)?-1:0;
    //     //cell is one behind
    //     cell_adjust = (sgn==1)?-1:0;

    //     ghost_cell[dir]+=gcell_adjust;
    //     cell_int[dir]+=cell_adjust;
    //     

    //     if(sgn == -1 && dir==1) 
    //     { // lo sides
    //         robin_a(ghost_cell) = 1.0;
    //         robin_b(ghost_cell) = 0.0;
    //         robin_f(ghost_cell) = prob_parm.V1/prob_parm.phis;
    //         bc_arr(ghost_cell) = prob_parm.V1/prob_parm.phis;
    //     }
    //     else if(sgn == 1 && dir==1)
    //     {
    //         robin_a(ghost_cell) = 1.0;
    //         robin_b(ghost_cell) = 0.0;
    //         robin_f(ghost_cell) = prob_parm.V2/prob_parm.phis;
    //         bc_arr(ghost_cell) = prob_parm.V2/prob_parm.phis;
    //     }
    //     else
    //     {
    //         robin_a(ghost_cell) = 0.0;
    //         robin_b(ghost_cell) = 1.0;
    //         robin_f(ghost_cell) = 0.0;
    //         bc_arr(ghost_cell) = 0.0;
    //     } 
    // }

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
                    amrex::Real Tg,amrex::Real Pg)

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
        
        //for all species
        if(sgn == -1) 
        { // lo sides
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
        }
        else
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
        } 

    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_vel(int i, int j, int k, int dir,
                     int specid,
                     Array4<Real> const& phi,
                     Array4<Real> const& efield,
                     Array4<Real> const& vel,
                     GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                     GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                     const GpuArray<int, AMREX_SPACEDIM> domlo,
                     const GpuArray<int, AMREX_SPACEDIM> domhi,
                     GpuArray<Real, AMREX_SPACEDIM> dx,
                     const Real time,
                     ProbParm const& prob_parm,
                     amrex::Real Tg, amrex::Real Pg)
    {
        IntVect leftcell(i, j, k);
        IntVect rightcell(i, j, k);
        IntVect face(i, j, k);
        leftcell[dir] -= 1;

        vel(i,j,k)=0.0;

        if(specid==E_ID)
        {
            //to impose 0 flux condition
            if( (face[dir]==domlo[dir]) || (face[dir]==(domhi[dir]+1)) )
            {
                vel(i,j,k)=0.0;
            }
            else
            {
                /*amrex::Real mu_left = mobility(left[0],left[1],left[2],specid,phi,prob_lo,
                  prob_hi,dx,time,prob_parm,Tg,Pg);

                  amrex::Real mu_right = mobility(right[0],right[1],right[2],specid,phi,prob_lo,
                  prob_hi,dx,time,prob_parm,Tg,Pg);

                  amrex::Real mu_face = 0.5*(mu_left+mu_right);*/

                Real efield_magL = std::sqrt(std::pow(phi(leftcell,EFX_ID),2.0) +
                                             std::pow(phi(leftcell,EFY_ID),2.0) +
                                             std::pow(phi(leftcell,EFZ_ID),2.0) );

                Real efield_magR = std::sqrt(std::pow(phi(rightcell,EFX_ID),2.0) +
                                             std::pow(phi(rightcell,EFY_ID),2.0) +
                                             std::pow(phi(rightcell,EFZ_ID),2.0) );

                Real efield_mag=0.5*(efield_magL+efield_magR)*prob_parm.phis;
                Real mu_face = -2.3987*std::pow(efield_mag,-0.26);

                Real efield_face = -(phi(rightcell,POT_ID) - phi(leftcell,POT_ID))/dx[dir]*prob_parm.phis;
                vel(i,j,k) = mu_face*efield_face;
            }
        }
    }
}
#endif
