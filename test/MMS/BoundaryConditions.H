#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>
#include <VarDefines.H>

using namespace amrex;
namespace plasmachem_transport
{
    AMREX_GPU_DEVICE AMREX_INLINE
    void potential_bc(int i, int j, int k,
                           int dir, int sgn,
                           Array4<Real> const& phi,
                           Array4<Real> const& bcarr,
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
        //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int(i,j,k);
        IntVect ghost_cell(i,j,k);
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        
        int gcell_adjust,cell_adjust;

        //ghost_cell is one behind
        gcell_adjust= (sgn==-1)?-1:0;
        //cell is one behind
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;
        
        if(sgn == -1) 
        { // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
        }
        else
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;
        } 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn, 
                    int spec_id, Array4<Real> const &phi, 
                    Array4<Real> const& bcarr,
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
        //note: bdryLo/bdryHi grabs the face indices from bx that are at the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int(i,j,k);
        IntVect ghost_cell(i,j,k);
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;

        int gcell_adjust,cell_adjust;

        gcell_adjust= (sgn==-1)?-1:0;
        cell_adjust = (sgn ==1)?-1:0;

        ghost_cell[dir]+=gcell_adjust;
        cell_int[dir]+=cell_adjust;
        
        if(sgn == -1) 
        { 
            // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 0.0;

        }
        else
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 1.0;
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
        Real xlen,ylen,zlen,maxlen;
        int bigdir=0;

        xlen=prob_hi[0]-prob_lo[0];
        ylen=prob_hi[1]-prob_lo[1];
        zlen=prob_hi[2]-prob_lo[2];


        IntVect face(i, j, k);
        IntVect leftcell(i, j, k);
        IntVect rightcell(i, j, k);
        leftcell[dir] -= 1;

        maxlen=std::max(xlen,std::max(ylen,zlen));

        if(maxlen==xlen) bigdir=0;
        if(maxlen==ylen) bigdir=1;
        if(maxlen==zlen) bigdir=2;

        vel(i,j,k)=0.0;
        
        if(dir==bigdir && specid==AR_ID)
        {
            //use this condition if using a robin bc that involves
            //velocity
            /*if( (right[dir]==domlo[dir]) || (right[dir]==(domhi[dir]+1)) )
              {
              vel(i,j,k,S1_ID)=0.0;
              }*/
            //else
            Real efield_magL = std::sqrt(std::pow(phi(leftcell,EFX_ID),2.0) +
                                         std::pow(phi(leftcell,EFY_ID),2.0) +
                                         std::pow(phi(leftcell,EFZ_ID),2.0) );
            
            Real efield_magR = std::sqrt(std::pow(phi(rightcell,EFX_ID),2.0) +
                                         std::pow(phi(rightcell,EFY_ID),2.0) +
                                         std::pow(phi(rightcell,EFZ_ID),2.0) );
            
            amrex::Real ndens = 0.0;
            amrex::Real mu_left = specMob(specid,phi(leftcell,ETEMP_ID),
                                           ndens, efield_magL,
                                           Tg);

            amrex::Real mu_right = specMob(specid,phi(rightcell,ETEMP_ID),
                                            ndens, efield_magR,
                                            Tg);

            amrex::Real mu_face = 0.5*(mu_left+mu_right);
            //amrex::Real Efld_face = 0.5*(phi(leftcell,EFX_ID+dir)+phi(rightcell,EFX_ID+dir));
            amrex::Real Efld_face = efield(face);
            vel(i,j,k) = mu_face*Efld_face;
        }
    }
}
#endif
