#include <Prob.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <PlasmaChem.H>
#include <UnivConstants.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    amrex::Vector<int> spec_chrg(NUM_SPECIES);

    void init()
    {
        // Get vector with species names
        CKSYMS_STR(specnames);

        // Get species charge numbers
        CKCHRG(&spec_chrg[0]);
    }    

    void close(){}
    
    int find_id(std::string specname){
      int loc=-1;
      auto it=std::find(specnames.begin(),specnames.end(),specname);
      if(it != specnames.end())
      {
          loc=it-specnames.begin();
      }
      return(loc);
    }
    
    int get_charge(int specid){
      return spec_chrg[specid];
    }
    amrex::Real get_molwt(int specid){
      return(1.0);
    }
}

namespace plasmachem_transport
{
    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real mobility(int specid,
                         Real Te, Real efield_x,
                         Real efield_y,Real efield_z,
                         GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                         GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                         GpuArray<Real, AMREX_SPACEDIM> dx,
                         const Real time,
                         ProbParm const& prob_parm,
                         Real Tg, Real Pg){
      // if(constant_transport){
      //
      // } else {
      //
      // }
      return(1.0);
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real diffusion_coeff(int specid,
                                Real Te, Real efield_x,
                                Real efield_y,Real efield_z,
                                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                GpuArray<Real, AMREX_SPACEDIM> dx,
                                const Real time,
                                ProbParm const& prob_parm,
                                Real Tg, Real Pg){
      // if(constant_transport){
      //
      // } else {
      //
      // }
      return(1.0);
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real collision_freq(int i, int j, int k,
                               int specid,
                               Array4<Real> const& phi,
                               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                               GpuArray<Real, AMREX_SPACEDIM> dx,
                               const Real time,
                               ProbParm const& prob_parm,
                               Real Tg, Real Pg){
      return(1.0);
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
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
                      amrex::Real Tg,amrex::Real Pg)
    {
      // TODO: Make this routine more flexible/general
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
            robin_f(ghost_cell) = prob_parm.V1;
            bc_arr(ghost_cell) = prob_parm.V1;
        }
        else if(sgn == 1)
        {
            robin_a(ghost_cell) = 1.0;
            robin_b(ghost_cell) = 0.0;
            robin_f(ghost_cell) = prob_parm.V2;
            bc_arr(ghost_cell) = prob_parm.V2;
        }
        else
        {
            robin_a(ghost_cell) = 0.0;
            robin_b(ghost_cell) = 1.0;
            robin_f(ghost_cell) = 0.0;
            bc_arr(ghost_cell) = 0.0;
        }
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
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
                    amrex::Real Tg,amrex::Real Pg){
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
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
                     amrex::Real Tg, amrex::Real Pg){
    }

}

namespace plasmachem_reactions
{
    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_react_source(int i, int j, int k,
                              int specid,
                              Array4<Real> const& phi,
                              GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                              GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                              GpuArray<Real, AMREX_SPACEDIM> dx,
                              const Real time,
                              ProbParm const& prob_parm,
                              amrex::Real Tg,amrex::Real Pg){

          // // Create array with species concentrations (1/m3 -> mol/cm3)
          // amrex::Real spec_C[NUM_SPECIES];
          // amrex::Real spec_wdot[NUM_SPECIES];
          // for(int sp=0; sp<NUM_SPECIES; sp++) spec_C[sp] = phi(i,j,k,sp) * 1.0e-6 / N_A;

          // // Get molar production rates
          // CKWC(Tg, spec_C, spec_wdot);

          // // Convert back to 1/m3
          // for()

          return(1.0);
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
    void compute_potential_source(int i, int j, int k,
                                  Array4<Real> const& phi,
                                  Array4<Real> const& source,
                                  GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                  GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                  GpuArray<Real, AMREX_SPACEDIM> dx,
                                  const Real time,
                                  ProbParm const& prob_parm){
    }

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_electron_inelastic_heating(int i, int j, int k,
                     Array4<Real> const& phi,
                     GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                     GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                     GpuArray<Real, AMREX_SPACEDIM> dx,
                     const Real time,
                     ProbParm const& prob_parm,
                     amrex::Real Tg,
                     amrex::Real pres){
      return(1.0);
    }

}
