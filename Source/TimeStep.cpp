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
#include <Transport.H>

// a wrapper for EstTimeStep
void Vidyut::ComputeDt()
{
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        dt[lev] = dt[lev - 1];
    }
}

// compute dt from CFL considerations
Real Vidyut::EstTimeStep(int lev)
{
    BL_PROFILE("Vidyut::EstTimeStep()");
    amrex::Real dt_est=fixed_dt;
    return dt_est;
}
