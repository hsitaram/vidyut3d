#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <Kernels_3d.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <Transport.H>
#include <Reactions.H>
#include <ProbParm.H>
#include <AMReX_MLABecLaplacian.H>

// advance solution to final time
void echemAMR::Evolve()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step + 1 << " starts ..." << std::endl;
    
        ComputeDt();

        if (max_level > 0 && regrid_int > 0)  // We may need to regrid
        {
            if (istep[0] % regrid_int == 0)
            {
                regrid(0, cur_time);
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << "ADVANCE with time = " << t_new[lev]
            << " dt = " << dt[0] << std::endl;
        }

        Vector< Array<MultiFab,AMREX_SPACEDIM> > flux(finest_level+1);
        Vector<MultiFab *> expl_src(finest_level+1);
        Vector<MultiFab *> Sborder(finest_level+1);
        
        //copy new to old and update time
        for(int lev=0;lev<=finest_level;lev++)
        {
            amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            t_old[lev] = t_new[lev];
            t_new[lev] += dt[0];
        }

        //allocate flux, expl_src, Sborder
        for(int lev=0;lev<=finest_level;lev++)
        {
            int num_grow=2;
            Sborder[lev] = new MultiFab(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder[lev]->setVal(0.0);
            FillPatch(lev, cur_time, *Sborder[lev], 0, Sborder[lev]->nComp());

            for (int lev = 0; lev <= finest_level; lev++)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    BoxArray ba = grids[lev];
                    ba.surroundingNodes(idim);
                    flux[lev][idim].define(ba, dmap[lev], 1, 0);
                    flux[lev][idim].setVal(0.0);
                }
                expl_src[lev]=new MultiFab(grids[lev], dmap[lev], 1, 0);
                expl_src[lev]->setVal(0.0);
            }
        }

        solve_potential(cur_time, Sborder);

        update_explsrc_at_all_levels(EDN_ID, Sborder, flux, expl_src, cur_time);
        implicit_solve_species(cur_time,dt[0],EDN_ID,Sborder,expl_src);

        if(elecenergy_solve)
        {
            update_explsrc_at_all_levels(EEN_ID, Sborder, flux, expl_src, cur_time);
            implicit_solve_species(cur_time,dt[0],EEN_ID, Sborder, expl_src);
        }

        for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
        {
            update_explsrc_at_all_levels(ind, Sborder, flux, expl_src, cur_time);
            implicit_solve_species(cur_time, dt[0], ind, Sborder, expl_src);
        }

        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
            ++istep[lev];

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }

        cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt[0] << std::endl;

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step + 1) % plot_int == 0)
        {
            last_plot_file_step = step + 1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step + 1) % chk_int == 0)
            {
                WriteCheckpointFile();
            }

            if (cur_time >= stop_time - 1.e-6 * dt[0]) break;
        }

        if (plot_int > 0 && istep[0] > last_plot_file_step)
        {
            WritePlotFile();
        }
    }
