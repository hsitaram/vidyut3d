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
        //mobility is scaled by charge
        amrex::Real mob=0.0;
        amrex::Real elecmob=-20.0;
        if(specid==EDN_ID)
        {
            mob=elecmob;
        }
        if(specid==I_ID)
        {
            mob=0.2;
        }
        if(specid==EEN_ID)
        {
            mob=fivebythree*elecmob;
        }

        return(mob);
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
        amrex::Real elecdcoeff=100.0;
        if(specid==EDN_ID)
        {
            dcoeff=elecdcoeff;
        }
        if(specid==I_ID)
        {
            dcoeff=0.01;
        }
        if(specid==EEN_ID)
        {
            dcoeff=fivebythree*elecdcoeff;
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
        amrex::Real mu = mobility(i,j,k,specid,phi,prob_lo,
                                   prob_hi,dx,time,prob_parm,Tg,Pg);

        amrex::Real mass = plasmachem::get_molwt(specid); 
        amrex::Real Z = plasmachem::get_charge(specid); 

        amrex::Real collfreq = Z*ECHARGE/mass/mu;

        return(collfreq);
    }
    
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
            robin_f(ghost_cell) = prob_parm.V1/prob_parm.phis;
            bc_arr(ghost_cell) = prob_parm.V1/prob_parm.phis;
        }
        else if(sgn == 1)
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = prob_parm.V2/prob_parm.phis;
            bc_arr(ghost_cell) = prob_parm.V2/prob_parm.phis;
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
        const int im1 = (dir == 0) ? i-1 : i;
        const int jm1 = (dir == 1) ? j-1 : j;
        const int km1 = (dir == 2) ? k-1 : k;

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
        
        //default
        robin_a(ghost_cell) = 0.0;
        robin_b(ghost_cell) = 1.0;
        robin_f(ghost_cell) = 0.0;

        if(spec_id==EDN_ID)
        {
            amrex::Real mob = mobility(cell_int[0],cell_int[1],cell_int[2],spec_id,phi,prob_lo,
                                       prob_hi,dx,time,prob_parm,Tg,Pg);
            amrex::Real Efield=phi(cell_int,EFX_ID+dir);
            amrex::Real driftvel=mob*Efield;

            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real D=diffusion_coeff(cell_int[0],cell_int[1],cell_int[2], 
                                          spec_id, phi, 
                                          prob_lo,prob_hi, dx, time, prob_parm,
                                          Tg,Pg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);

            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);

            robin_a(ghost_cell) = driftvel;
            robin_b(ghost_cell) = -D;
            robin_f(ghost_cell) = ne*cbar/4.0;
        }

        if(spec_id==I_ID)
        {
            amrex::Real mob = mobility(cell_int[0],cell_int[1],cell_int[2],spec_id,phi,prob_lo,
                                       prob_hi,dx,time,prob_parm,Tg,Pg);
            amrex::Real Efield=phi(cell_int,EFX_ID+dir);

            amrex::Real driftvel=mob*Efield;

            //default setting
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;

            //directed at the surface
            if(driftvel*outward_normal[dir] >= 0.0)
            {
                /*amrex::Prit()<<"ion neumann,sgn:"<<sgn<<"\t"<<driftvel<<"\t"<<
                  mob<<"\t"<<Efield<<"\n";*/
                //drift dominated flux and quenching
                //mu n E - D grad n = mu n E
                //implied D grad n = 0
                robin_a(ghost_cell) = 0.0;
                robin_b(ghost_cell) = 1.0;
                robin_f(ghost_cell) = 0.0;
            }
            else
            {
                //amrex::Print()<<"ion robin,sgn:"<<sgn<<"\n";
                amrex::Real D=diffusion_coeff(cell_int[0],cell_int[1],cell_int[2], 
                                              spec_id, phi, 
                                              prob_lo,prob_hi, dx, time, prob_parm,
                                              Tg,Pg);
                //no flux
                //mu n E - D grad n = 0
                robin_a(ghost_cell) = driftvel;
                robin_b(ghost_cell) = -D;
                robin_f(ghost_cell) = 0.0;
                /*robin_a(ghost_cell) = 0.0;
                  robin_b(ghost_cell) = 1.0;
                  robin_f(ghost_cell) = 0.0;*/
            }
        }

        if(spec_id==EEN_ID)
        {
            amrex::Real mob = mobility(cell_int[0],cell_int[1],cell_int[2],spec_id,phi,prob_lo,
                                       prob_hi,dx,time,prob_parm,Tg,Pg);
            amrex::Real Efield=phi(cell_int,EFX_ID+dir);
            amrex::Real driftvel=mob*Efield;

            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real D=diffusion_coeff(cell_int[0],cell_int[1],cell_int[2], 
                                          spec_id, phi, 
                                          prob_lo,prob_hi, dx, time, prob_parm,
                                          Tg,Pg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);

            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);

            robin_a(ghost_cell) = driftvel;
            robin_b(ghost_cell) = -D;
            robin_f(ghost_cell) = 2.0*K_B*Te*ne*cbar/4.0;
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
                     amrex::Real Tg,amrex::Real Pg)
    {
        IntVect left(i, j, k);
        IntVect right(i, j, k);
        left[dir] -= 1;

        vel(i,j,k)=0.0;
        amrex::Real mu_left = mobility(left[0],left[1],left[2],specid,phi,prob_lo,
                                       prob_hi,dx,time,prob_parm,Tg,Pg);

        amrex::Real mu_right = mobility(right[0],right[1],right[2],specid,phi,prob_lo,
                                        prob_hi,dx,time,prob_parm,Tg,Pg);

        amrex::Real mu_face = 0.5*(mu_left+mu_right);

        Real efield_face = -(phi(right,POT_ID) - phi(left,POT_ID))/dx[dir];

        if( (right[dir]==domlo[dir]) || (right[dir]==(domhi[dir]+1)) )
        {
            vel(i,j,k)=0.0;
        }
        else
        {
            vel(i,j,k) = mu_face*efield_face;
        }
    }


}
#endif
