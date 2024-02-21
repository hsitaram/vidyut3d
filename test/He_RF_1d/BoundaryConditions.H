#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <Chemistry.H>

using namespace amrex;
namespace user_transport
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
                      Real Tg, Real Pg, Real app_voltage)
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
            robin_f(ghost_cell) = app_voltage;
            bc_arr(ghost_cell)  = app_voltage;
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
        const int im1 = (dir == 0) ? i-1 : i;
        const int jm1 = (dir == 1) ? j-1 : j;
        const int km1 = (dir == 2) ? k-1 : k;

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

        if(spec_id==E_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            amrex::Real ndens = 0.0;
            for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cell_int,sp);
            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);
            
            amrex::Real Te=phi(cell_int,ETEMP_ID);
            amrex::Real cbar=sqrt(8.0*K_B*Te/PI/ME);
            amrex::Real flx = ne*cbar/4.0*Real(sgn);

            bc_arr(ghost_cell)=(-flx/D);
        }

        if(spec_id==HEp_ID || spec_id==HE2p_ID)
        {
            bc_arr(ghost_cell)=0.0;
        }

        if(spec_id==HEm_ID || spec_id==HE2m_ID)
        {
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            
            amrex::Real ndens = 0.0;
            for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cell_int,sp);
            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);
            
            // amrex::Real cbar=sqrt(8.0*K_B*Tg/PI/plasmachem::get_molwt(spec_id));
            amrex::Real cbar=sqrt(8.0*K_B*Tg/PI/(h_global_mw[spec_id]*M_AMU));
            amrex::Real flx = phi(i,j,k,spec_id)*cbar/4.0*Real(sgn);
            bc_arr(ghost_cell)=(-flx/D);
        }

        if(spec_id==EEN_ID)
        {
            amrex::Real ne=phi(cell_int,E_ID);
            amrex::Real efieldmag=std::sqrt(std::pow(phi(cell_int,EFX_ID),2.0)+
                                            std::pow(phi(cell_int,EFY_ID),2.0)+
                                            std::pow(phi(cell_int,EFZ_ID),2.0));
            
            amrex::Real Te=phi(cell_int,ETEMP_ID);
            
            amrex::Real ndens = 0.0;
            for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(cell_int,sp);
            amrex::Real D=specDiff(spec_id, phi(cell_int,ETEMP_ID), ndens, efieldmag, Tg);
            
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
        amrex::Real efieldmag_L=std::sqrt(std::pow(phi(leftcell,EFX_ID),2.0)+
                                        std::pow(phi(leftcell,EFY_ID),2.0)+
                                        std::pow(phi(leftcell,EFZ_ID),2.0));
        
        amrex::Real efieldmag_R=std::sqrt(std::pow(phi(rightcell,EFX_ID),2.0)+
                                        std::pow(phi(rightcell,EFY_ID),2.0)+
                                        std::pow(phi(rightcell,EFZ_ID),2.0));
        
        amrex::Real ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(leftcell,sp);
        amrex::Real mu_leftcell = specMob(specid,phi(leftcell,ETEMP_ID), ndens, efieldmag_L, Tg);

        ndens = 0.0;
        for(int sp=0; sp<NUM_SPECIES; sp++) ndens += phi(rightcell,sp);
        amrex::Real mu_rightcell = specMob(specid,phi(rightcell,ETEMP_ID), ndens, efieldmag_R, Tg);

        amrex::Real mu_face = 0.5*(mu_leftcell+mu_rightcell);

        //Real efield_face = -(phi(rightcell,POT_ID) - phi(leftcell,POT_ID))/dx[dir];
        //Real efield_face = 0.5*(phi(rightcell,EFX_ID+dir) + phi(leftcell,EFX_ID+dir));
        Real efield_face = efield(face);
        Real driftvel = mu_face*efield_face;

        if(dir==dischargedir)
        {
            if( (face[dir]==domlo[dir]) || (face[dir]==(domhi[dir]+1)) )
            {
                Real outward_normal_sgn=(face[dir]==domlo[dir])?-1.0:1.0;

                if( (specid==E_ID) || (specid==EEN_ID) )
                {
                    //false Neumann for robin bc
                    vel(face)=0.0;
                }
                else if(specid==HEp_ID || specid==HE2p_ID)
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
