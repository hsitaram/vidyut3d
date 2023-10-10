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
    amrex::Real mobility(int i, int j, int k,
                         int specid,
                         Array4<Real> const& phi,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         const Real time,
                         ProbParm const& prob_parm,
                         Real Tg, Real Pg)

    {
        return(1.0);
    }
    
    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real diffusion_coeff(int i, int j, int k,
                                int specid,
                                Array4<Real> const& phi,
                                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                GpuArray<Real, AMREX_SPACEDIM> dx,
                                const Real time,
                                ProbParm const& prob_parm,
                                Real Tg, Real Pg)

    {
        amrex::Real dcoeff=0.0;
        if(specid==S1_ID)
        {
           dcoeff = 0.1;
        }
        else
        {
          dcoeff = 1.0;
        }
        return(dcoeff);
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real collision_freq(int i, int j, int k,
                               int specid,
                               Array4<Real> const& phi,
                               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                               GpuArray<Real, AMREX_SPACEDIM> dx,
                               const Real time,
                               ProbParm const& prob_parm,
                               Real Tg, Real Pg)

    {
        return(1.0);
    }
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
        //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
        //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
        //so the ghost cell index at left side is i-1 while it is i on the right
        IntVect cell_int(i,j,k);
        IntVect ghost_cell(i,j,k);
        amrex::Real outward_normal[AMREX_SPACEDIM]={0.0};
        outward_normal[dir]=sgn;
        
        //low side
        if(sgn==-1)
        {
           //cell is i,j,k
           //ghost_cell is one behind
           ghost_cell[dir]-=1;
        }
        
        //high side
        if(sgn==1)
        {
           //ghost_cell is i,j,k
           //cell is one behind
           cell_int[dir]-=1;
        }

        if(sgn == -1) 
        { // lo sides
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = 1.0;
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
        
        //low side
        if(sgn==-1)
        {
           //cell is i,j,k
           //ghost_cell is one behind
           ghost_cell[dir]-=1;
        }
        
        //high side
        if(sgn==1)
        {
           //ghost_cell is i,j,k
           //cell is one behind
           cell_int[dir]-=1;
        }
        
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
            if(spec_id==S1_ID)
            {
                robin_f(ghost_cell) = 1.0;
            }
            else //also for elec energy and density
            {
                robin_f(ghost_cell) = 0.0;
            }
        } 

    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_vel(int i, int j, int k, int dir,
                     int specid,
                     Array4<Real> const& phi,
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


        IntVect left(i, j, k);
        IntVect right(i, j, k);
        left[dir] -= 1;

        maxlen=std::max(xlen,std::max(ylen,zlen));

        if(maxlen==xlen) bigdir=0;
        if(maxlen==ylen) bigdir=1;
        if(maxlen==zlen) bigdir=2;

        vel(i,j,k)=0.0;

        if(bigdir==dir && specid==S1_ID)
        {
            vel(i,j,k)
            =0.5*(phi(left,EFX_ID+dir)+phi(right,EFX_ID+dir));
        }
    }
}
#endif
