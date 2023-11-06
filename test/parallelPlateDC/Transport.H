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
            //mob=0.0;
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
            //dcoeff=0.0;
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
            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real D=diffusion_coeff(cell_int[0],cell_int[1],cell_int[2], 
                                          spec_id, phi, 
                                          prob_lo,prob_hi, dx, time, prob_parm,
                                          Tg,Pg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            amrex::Real flx = ne*cbar/4.0*Real(sgn);

            bc_arr(ghost_cell)=(-flx/D);
        }

        if(spec_id==I_ID)
        {
            bc_arr(ghost_cell)=0.0;
        }

        if(spec_id==EEN_ID)
        {
            amrex::Real ne=phi(cell_int,EDN_ID);
            amrex::Real D=diffusion_coeff(cell_int[0],cell_int[1],cell_int[2], 
                                          spec_id, phi, 
                                          prob_lo,prob_hi, dx, time, prob_parm,
                                          Tg,Pg);

            amrex::Real Te=phi(cell_int,ETEMP_ID);

            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            
            amrex::Real flx = (2.0*K_B*Te)*ne*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
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
                     amrex::Real Tg,amrex::Real Pg)
    {
        IntVect leftcell(i, j, k);
        IntVect rightcell(i, j, k);
        IntVect face(i, j, k);
        leftcell[dir] -= 1;
        
        //in this case, discharge is along x
        int dischargedir=0;

        vel(face)=0.0;
        amrex::Real mu_leftcell = mobility(leftcell[0],leftcell[1],leftcell[2],specid,phi,prob_lo,
                                       prob_hi,dx,time,prob_parm,Tg,Pg);

        amrex::Real mu_rightcell = mobility(rightcell[0],rightcell[1],rightcell[2],specid,phi,prob_lo,
                                        prob_hi,dx,time,prob_parm,Tg,Pg);

        amrex::Real mu_face = 0.5*(mu_leftcell+mu_rightcell);

        //Real efield_face = -(phi(rightcell,POT_ID) - phi(leftcell,POT_ID))/dx[dir];
        Real efield_face = 0.5*(phi(rightcell,EFX_ID+dir) + phi(leftcell,EFX_ID+dir));
        Real driftvel = mu_face*efield_face;

        if(dir==dischargedir)
        {
            if( (face[dir]==domlo[dir]) || (face[dir]==(domhi[dir]+1)) )
            {
                Real outward_normal_sgn=(face[dir]==domlo[dir])?-1.0:1.0;

                if( (specid==EDN_ID) || (specid==EEN_ID) )
                {
                    //false Neumann for robin bc
                    vel(face)=0.0;
                }
                else if(specid==I_ID)
                {
                    //if driftvel is pointing towards the wall
                    //then flux=v n - D dndx = v n => dndx=0 true Neumann bc
                    //else  flux = v n - D dndx = 0 false Neumann for robin bc
                    //
                    //Real iondriftvel = mu_face*0.5*(phi(rightcell,EFX_ID+dir) + phi(leftcell,EFX_ID+dir));
                    Real iondriftvel=driftvel;
                    vel(face)=(iondriftvel*outward_normal_sgn > 0)?iondriftvel:0.0;
                    
                    //vel(face)=0.0;
                    //amrex::Print()<<"ion vel,sgn,face:"<<vel(face)<<"\t"<<outward_normal_sgn<<"\t"<<face<<"\n";
                }
            }
            else
            {
                vel(face) = driftvel;
            }
        }
        else
        {
            vel(face) = driftvel;
        }
    }
}
#endif
