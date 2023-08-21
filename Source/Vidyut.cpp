#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <Kernels_3d.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <ProbParm.H>
#include <stdio.h>
#include <Constants.H>

using namespace amrex;

ProbParm* echemAMR::h_prob_parm = nullptr;
ProbParm* echemAMR::d_prob_parm = nullptr;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
echemAMR::echemAMR()
{
    ReadParameters();
    h_prob_parm = new ProbParm{};
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));
    amrex_probinit(*h_prob_parm, *d_prob_parm);

    plasma_param_names.resize(NUM_PLASMAVARS);
    plasma_param_names[0]="Electron_density";
    plasma_param_names[1]="Electron_energy";
    plasma_param_names[2]="Eden_gradx";
    plasma_param_names[3]="Eden_grady";
    plasma_param_names[4]="Eden_gradz";
    plasma_param_names[5]="Potential";
    plasma_param_names[6]="Efieldx";
    plasma_param_names[7]="Efieldy";
    plasma_param_names[8]="Efieldz";
    
    allvarnames.resize(NVAR);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
        allvarnames[i] = plasmachem::specnames[i];
    }
    for(int i=0;i<NUM_PLASMAVARS;i++)
    {
        allvarnames[i+NUM_SPECIES] = plasma_param_names[i];
    }

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev)
    {
        nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    for(int lev=0;lev<nlevs_max;lev++)
    {
        dt[lev]=fixed_dt;
    }

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    
    ParmParse pp("vidyut");
    pp.queryarr("pot_bc_lo", pot_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("pot_bc_hi", pot_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("eden_bc_lo", eden_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("eden_bc_hi", eden_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("eenrg_bc_lo", eenrg_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("eenrg_bc_hi", eenrg_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("ion_bc_lo", ion_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("ion_bc_hi", ion_bc_hi, 0, AMREX_SPACEDIM);
    
    pp.queryarr("neutral_bc_lo", neutral_bc_lo, 0, AMREX_SPACEDIM);
    pp.queryarr("neutral_bc_hi", neutral_bc_hi, 0, AMREX_SPACEDIM);

    bcspec.resize(NVAR);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {

        int bctype=(geom[0].isPeriodic(idim))?BCType::int_dir:BCType::foextrap;

        for (int sp=0; sp < NVAR; sp++) 
        {
            bcspec[sp].setLo(idim, bctype);
            bcspec[sp].setHi(idim, bctype);
        }
    }

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
}

echemAMR::~echemAMR()
{
    delete h_prob_parm;
    The_Arena()->free(d_prob_parm);
}
// initializes multilevel data
void echemAMR::InitData()
{
    ProbParm* localprobparm = d_prob_parm;

    if (restart_chkfile == "")
    {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0)
        {
            WriteCheckpointFile();
        }

    } 
    else
    {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0)
    {
        WritePlotFile();
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void echemAMR::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;

    // only do this during the first call to ErrorEst
    if (first)
    {
        first = false;
        // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
        ParmParse pp("vidyut");
            if (pp.contains("tagged_vars"))
        {
            int nvars = pp.countval("tagged_vars");
            refine_phi.resize(nvars);
            refine_phigrad.resize(nvars);
            refine_phi_comps.resize(nvars);
            std::string varname;
            for (int i = 0; i < nvars; i++)
            {
                pp.get("tagged_vars", varname, i);
                pp.get((varname + "_refine").c_str(), refine_phi[i]);
                pp.get((varname + "_refinegrad").c_str(), refine_phigrad[i]);

                int spec_id = plasmachem::find_id(varname);
                if (spec_id == -1)
                {
                    int varname_id=-1;
                    auto it=std::find(plasma_param_names.begin(),plasma_param_names.end(),varname);
                    if(it != plasma_param_names.end())
                    {
                        varname_id=it-plasma_param_names.begin();
                    }

                    if (varname_id == -1)
                    {
                        Print() << "Variable name:" << varname << " not found for tagging\n";
                        amrex::Abort("Invalid tagging variable");
                    }
                    else
                    {
                        refine_phi_comps[i] = varname_id+NUM_SPECIES;
                    }
                }
                else
                {
                   refine_phi_comps[i] = spec_id;
                }
            }
        }
    }

    if (refine_phi.size() == 0) return;

    //    const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;

    const MultiFab& state = phi_new[lev];
    MultiFab Sborder(grids[lev], dmap[lev], state.nComp(), 1);
    FillPatch(lev, time, Sborder, 0, Sborder.nComp());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto statefab = Sborder.array(mfi);
            const auto tagfab = tags.array(mfi);

            amrex::Real* refine_phi_dat = refine_phi.data();
            amrex::Real* refine_phigrad_dat = refine_phigrad.data();
            int* refine_phi_comps_dat = refine_phi_comps.data();
            int ntagged_comps = refine_phi_comps.size();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                state_based_refinement(i, j, k, tagfab, statefab, refine_phi_dat, refine_phi_comps_dat, ntagged_comps, tagval);
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                stategrad_based_refinement(i, j, k, tagfab, statefab, refine_phigrad_dat, refine_phi_comps_dat, ntagged_comps, tagval);
            });
        }
    }
}

// read in some parameters from inputs file
void echemAMR::ReadParameters()
{
    {
        ParmParse pp; // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.
        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart", restart_chkfile);
    }

    {
        ParmParse pp("vidyut");

        pp.query("dt", fixed_dt);
        pp.query("use_hypre",use_hypre);

        pp.query("linsolve_reltol",linsolve_reltol);
        pp.query("linsolve_abstol",linsolve_abstol);
        pp.query("linsolve_bot_reltol",linsolve_bot_reltol);
        pp.query("linsolve_bot_abstol",linsolve_bot_abstol);

        pp.query("linsolve_num_pre_smooth",linsolve_num_pre_smooth);
        pp.query("linsolve_num_post_smooth",linsolve_num_post_smooth);
        pp.query("linsolve_num_final_smooth",linsolve_num_final_smooth);
        pp.query("linsolve_num_bottom_smooth",linsolve_num_bottom_smooth);

        pp.query("linsolve_maxiter",linsolve_maxiter);
        pp.query("linsolve_max_coarsening_level",linsolve_max_coarsening_level);
        pp.query("bound_specden", bound_specden);
        pp.query("min_species_density",min_species_density);
        pp.query("elecenergy_solve",elecenergy_solve);
    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void echemAMR::GetData(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    } else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    } else
    {
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}
