using namespace amrex;

namespace plasmachem
{
    extern amrex::Vector<std::string> specnames;
    extern amrex::Vector<int> spec_chrg;

    void init();
    void close();
    int find_id(std::string specname);
    
    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real get_bg_molwt(int i, int j, int k,
                     Array4<Real> const& phi,
                     ProbParm const& prob_parm)
    {
        amrex::Real molwt=40.0*M_AMU;    
        return(molwt);
    }
    
    AMREX_GPU_HOST_DEVICE int get_charge(int specid);
    AMREX_GPU_HOST_DEVICE amrex::Real get_molwt(int specid);
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
                         Real Tg, Real Pg);

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real diffusion_coeff(int specid,
                                Real Te, Real efield_x,
                                Real efield_y,Real efield_z,
                                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                GpuArray<Real, AMREX_SPACEDIM> dx,
                                const Real time,
                                ProbParm const& prob_parm,
                                Real Tg, Real Pg);

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real collision_freq(int i, int j, int k,
                               int specid,
                               Array4<Real> const& phi,
                               GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                               GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                               GpuArray<Real, AMREX_SPACEDIM> dx,
                               const Real time,
                               ProbParm const& prob_parm,
                               Real Tg, Real Pg);

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
                      amrex::Real Tg,amrex::Real Pg);

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
                    amrex::Real Tg,amrex::Real Pg);

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
                     amrex::Real Tg, amrex::Real Pg);
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
                              amrex::Real Tg,amrex::Real Pg);

    // AMREX_GPU_DEVICE AMREX_INLINE
    void compute_potential_source(int i, int j, int k,
                                  Array4<Real> const& phi,
                                  Array4<Real> const& source,
                                  GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                                  GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                                  GpuArray<Real, AMREX_SPACEDIM> dx,
                                  const Real time,
                                  ProbParm const& prob_parm);

    // AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real compute_electron_inelastic_heating(int i, int j, int k,
                     Array4<Real> const& phi,
                     GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                     GpuArray<Real, AMREX_SPACEDIM> prob_hi,
                     GpuArray<Real, AMREX_SPACEDIM> dx,
                     const Real time,
                     ProbParm const& prob_parm,
                     amrex::Real Tg,
                     amrex::Real pres);

}
