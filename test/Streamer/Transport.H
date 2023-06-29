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
                           Array4<Real> const& robin_a,
                           Array4<Real> const& robin_b,
                           Array4<Real> const& robin_f,
                           GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                           GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                           GpuArray<Real, AMREX_SPACEDIM> dx,
                           const Real time,
                           ProbParm const& prob_parm)
    {
        const int im1 = (dir == 0) ? i-1 : i;
        const int jm1 = (dir == 1) ? j-1 : j;
        const int km1 = (dir == 2) ? k-1 : k;
        
        if(sgn == -1) 
        { // lo sides
            robin_a(im1,jm1,km1) = 1.0;
            robin_b(im1,jm1,km1) = 0.0;
            robin_f(im1,jm1,km1) = prob_parm.V1;
        }
        else
        {
            robin_a(i,j,k) = 1.0;
            robin_b(i,j,k) = 0.0;
            robin_f(i,j,k) = prob_parm.V2;
        } 
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void species_bc(int i,int j,int k, int dir, int sgn, 
                             int spec_id, Array4<Real> const &phi, 
                              Array4<Real> const& robin_a,
                              Array4<Real> const& robin_b,
                             Array4<Real> const& robin_f,
                             GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                            GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                               const Real time,
                              ProbParm const& prob_parm)
                             
    {
        const int im1 = (dir == 0) ? i-1 : i;
        const int jm1 = (dir == 1) ? j-1 : j;
        const int km1 = (dir == 2) ? k-1 : k;
        
        if(sgn == -1) 
        { // lo sides
            robin_a(im1,jm1,km1) = 0.0;
            robin_b(im1,jm1,km1) = 1.0;
            robin_f(im1,jm1,km1) = 0.0;
        }
        else
        {
            robin_a(i,j,k) = 0.0;
            robin_b(i,j,k) = 1.0;
            robin_f(i,j,k) = 0.0;
        } 

    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_dcoeff(int i, int j, int k,
                        Array4<Real> const& phi,
                        Array4<Real> const& dcoeff,
                        GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                        GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                        GpuArray<Real, AMREX_SPACEDIM> dx,
                        const Real time,
                        ProbParm const& prob_parm)
    {
        Real elecfield_mag=std::sqrt(std::pow(phi(i,j,k,EFX_ID),2.0) +
                                     std::pow(phi(i,j,k,EFY_ID),2.0) +
                                     std::pow(phi(i,j,k,EFZ_ID),2.0) );
        
        dcoeff(i,j,k,E_ID)   = 4.3628e-3 * std::pow(elecfield_mag,0.22);
        dcoeff(i,j,k,I_ID)   = 0.0;
    }

    AMREX_GPU_DEVICE AMREX_INLINE
    void compute_vel(int i, int j, int k, int dir,
                      Array4<Real> const& phi,
                      Array4<Real> const& vel,
                      GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                      GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                      const GpuArray<int, AMREX_SPACEDIM> domlo,
                      const GpuArray<int, AMREX_SPACEDIM> domhi,
                      GpuArray<Real, AMREX_SPACEDIM> dx,
                      const Real time,
                      ProbParm const& prob_parm)
    {
        vel(i,j,k,I_ID) = 0.0;
        IntVect left(i, j, k);
        IntVect right(i, j, k);
        left[dir] -= 1;

        if( (right[dir]==domlo[dir]) || (right[dir]==(domhi[dir]+1)) )
        {
            vel(i,j,k,E_ID)=0.0;
        }
        else
        {

           Real efield_magL = std::sqrt(std::pow(phi(left,EFX_ID),2.0) +
                                        std::pow(phi(left,EFY_ID),2.0) +
                                        std::pow(phi(left,EFZ_ID),2.0) );
           
           Real efield_magR = std::sqrt(std::pow(phi(right,EFX_ID),2.0) +
                                        std::pow(phi(right,EFY_ID),2.0) +
                                        std::pow(phi(right,EFZ_ID),2.0) );

            Real efield_mag=0.5*(efield_magL+efield_magR);
           
            Real efield = -(phi(right,POT_ID) - phi(left,POT_ID))/dx[dir];
            Real mobility = -2.3987*std::pow(efield_mag,-0.26); 

            vel(i,j,k,E_ID) = mobility*efield;
        }
    }
}
#endif
