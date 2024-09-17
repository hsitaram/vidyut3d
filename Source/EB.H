#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#ifdef AMREX_USE_EB
#include <AMReX_EB2_IF_Cylinder.H>

#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

#include <AMReX_EB_utils.H>
#endif

namespace EBtools
{
    extern amrex::EBFArrayBoxFactory* ebfactory;
    extern amrex::MultiFab * lsphi;
    extern int ls_refinement;
    extern bool using_levelset_geometry;
    void init_eb(const amrex::Geometry &geom,const amrex::BoxArray &ba,const amrex::DistributionMapping &dm, int max_lev, int nghost);
}
