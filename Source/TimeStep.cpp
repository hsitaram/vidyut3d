#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Prob.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <BoundaryConditions.H>

// a wrapper for EstTimeStep
void Vidyut::ComputeDt()
{
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        dt[lev] = dt[lev - 1];
    }
}

// compute dt from CFL considerations
void Vidyut::find_time_scales(int lev,amrex::Real& dt_edrift,amrex::Real &dt_ediff,
            amrex::Real& dt_diel_relax)
{
    BL_PROFILE("Vidyut::EstTimeStep()");
    dt_edrift = std::numeric_limits<Real>::max();
    dt_diel_relax = std::numeric_limits<Real>::max();
    dt_ediff = std::numeric_limits<Real>::max();

    const auto dx = geom[lev].CellSizeArray();
    MultiFab& S_new = phi_new[lev];
    
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;

    MultiFab edriftvel(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab ediff(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
    MultiFab mue_ne(S_new.boxArray(), S_new.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> state_array = S_new.array(mfi);
            Array4<Real> edriftvel_array = edriftvel.array(mfi);
            Array4<Real> ediff_array = ediff.array(mfi);
            Array4<Real> mue_ne_array = mue_ne.array(mfi);
            auto prob_lo = geom[lev].ProbLoArray();
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                amrex::Real etemp=state_array(i,j,k,ETEMP_ID); 
                amrex::Real Esum = 0.0;
                for (int dim = 0; dim < AMREX_SPACEDIM; dim++) Esum += std::pow(state_array(i,j,k,EFX_ID+dim),2.0);
                amrex::Real efield_mag=std::sqrt(Esum);
               
                amrex::Real ndens = 0.0; 
                for(int sp=0; sp<NUM_SPECIES; sp++) ndens += state_array(i,j,k,sp);

                amrex::Real dcoeff = (const_ele_trans) ? ele_diff/ndens:specDiff(E_IDX, etemp, ndens,
                                               efield_mag,captured_gastemp);
                
                amrex::Real mu = (const_ele_trans) ? ele_mob/ndens:specMob(E_IDX, etemp, ndens,
                                               efield_mag, captured_gastemp);
                
                edriftvel_array(i,j,k)=amrex::Math::abs(mu)*efield_mag;
                ediff_array(i,j,k)=dcoeff;
                mue_ne_array(i,j,k)=amrex::Math::abs(mu)*state_array(i,j,k,E_IDX);
            });
        }
    }

    amrex::Real max_edriftvel = edriftvel.norm0(0,0,true);
    amrex::Real max_ediff    = ediff.norm0(0,0,true);
    amrex::Real max_muene    = mue_ne.norm0(0,0,true);
    if(max_edriftvel > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_edrift = std::min(dt_edrift, dx[i] / max_edriftvel);
        }
    }
    if(max_ediff > 0)
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) 
        {
            dt_ediff = std::min(dt_ediff, 0.5*dx[i]*dx[i]/max_ediff/AMREX_SPACEDIM);
        }
    }
    if(max_muene > 0)
    {
        dt_diel_relax=EPS0/ECHARGE/max_muene;
    }

    ParallelDescriptor::ReduceRealMin(dt_edrift);
    ParallelDescriptor::ReduceRealMin(dt_ediff);
    ParallelDescriptor::ReduceRealMin(dt_diel_relax);
}
