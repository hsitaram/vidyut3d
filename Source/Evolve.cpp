#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>
#include <ProbParm.H>
#include <Vidyut.H>
#include <Chemistry.H>
#include <PlasmaChem.H>
// #include <Transport.H>
// #include <Reactions.H>
#include <compute_flux_3d.H>
#include <AMReX_MLABecLaplacian.H>

// advance solution to final time
void Vidyut::Evolve()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    int plotfilenum=0;
    int chkfilenum=0;

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
        amrex::Real dt_common=dt[0]; //no subcycling
        int num_grow=2;

        Vector< Array<MultiFab,AMREX_SPACEDIM> > flux(finest_level+1);

        //face centered efield and electron density gradient
        Vector< Array<MultiFab,AMREX_SPACEDIM> > efield_fc(finest_level+1);
        Vector< Array<MultiFab,AMREX_SPACEDIM> > gradne_fc(finest_level+1);
        Vector< Array<MultiFab,AMREX_SPACEDIM> > grad_fc(finest_level+1);
        Vector<MultiFab> expl_src(finest_level+1);
        Vector<MultiFab> rxn_src(finest_level+1);
        Vector<MultiFab> Sborder(finest_level+1);

        //copy new to old and update time
        for(int lev=0;lev<=finest_level;lev++)
        {
            amrex::MultiFab::Copy(phi_old[lev], phi_new[lev], 
                                  0, 0, phi_new[lev].nComp(), 0);
            t_old[lev] = t_new[lev];
            t_new[lev] += dt_common;
        }

        //allocate flux, expl_src, Sborder
        for(int lev=0;lev<=finest_level;lev++)
        {
            Sborder[lev].define(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
            Sborder[lev].setVal(0.0);
            FillPatch(lev, cur_time, Sborder[lev], 0, Sborder[lev].nComp());

            for (int lev = 0; lev <= finest_level; lev++)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    BoxArray ba = grids[lev];
                    ba.surroundingNodes(idim);

                    flux[lev][idim].define(ba, dmap[lev], 1, 0);
                    flux[lev][idim].setVal(0.0);
                    
                    efield_fc[lev][idim].define(ba, dmap[lev], 1, 0);
                    efield_fc[lev][idim].setVal(0.0);
                    
                    gradne_fc[lev][idim].define(ba, dmap[lev], 1, 0);
                    gradne_fc[lev][idim].setVal(0.0);
                    
                    grad_fc[lev][idim].define(ba, dmap[lev], 1, 0);
                    grad_fc[lev][idim].setVal(0.0);
                }
                expl_src[lev].define(grids[lev], dmap[lev], 1, 0);
                expl_src[lev].setVal(0.0);

                rxn_src[lev].define(grids[lev], dmap[lev], NUM_SPECIES+1, 0);
                rxn_src[lev].setVal(0.0);
            }
        }

        solve_potential(cur_time, Sborder, pot_bc_lo, pot_bc_hi, efield_fc);
        
        // Calculate the reactive source terms for all species/levels
        if(do_reactions){
            update_rxnsrc_at_all_levels(Sborder, rxn_src, cur_time);
        }

        // note that phi_new is updated instead of sborder
        // so older potential and efield are used as opposed to new ones
        // call fillpatch to improve implicitness
        /*for(int lev=0;lev<=finest_level;lev++)
          {
          FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
          }*/

        for(unsigned int ind=0;ind<NUM_SPECIES;ind++)
        {
            update_explsrc_at_all_levels(ind, Sborder, flux, rxn_src, efield_fc, expl_src, cur_time);

            //electrons and ions
            if(plasmachem::get_charge(ind)!=0)
            {
                if(ind == E_IDX){
                    implicit_solve_scalar(cur_time,dt_common,E_IDX,Sborder,expl_src,eden_bc_lo,eden_bc_hi, gradne_fc);
                } else {
                    implicit_solve_scalar(cur_time, dt_common, ind, Sborder, expl_src,ion_bc_lo,ion_bc_hi, grad_fc);
                }
            }
            //neutrals
            else
            {
                implicit_solve_scalar(cur_time, dt_common, ind, Sborder, expl_src,neutral_bc_lo,neutral_bc_hi, grad_fc);
            }

            /*for(int lev=0;lev<=finest_level;lev++)
              {
              FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
              }*/
        }

        /*for(int lev=0;lev<=finest_level;lev++)
          {
          FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
          }*/

        if(elecenergy_solve)
        {
            update_explsrc_at_all_levels(EEN_ID, Sborder, flux, rxn_src, efield_fc, expl_src, cur_time);
            for (int lev = 0; lev <= finest_level; lev++)
            {
                compute_elecenergy_source(lev, num_grow, Sborder[lev], 
                                          efield_fc[lev], gradne_fc[lev], rxn_src[lev],
                                          expl_src[lev], cur_time, dt_common);
            }
            implicit_solve_scalar(cur_time,dt_common,EEN_ID, Sborder, 
                                  expl_src,eenrg_bc_lo,eenrg_bc_hi, grad_fc);

            /*for(int lev=0;lev<=finest_level;lev++)
              {
              FillPatch(lev, cur_time+dt_common, Sborder[lev], 0, Sborder[lev].nComp());
              }*/
        }

        AverageDown ();

        for (int lev = 0; lev <= finest_level; lev++)
            ++istep[lev];

        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }

        cur_time += dt_common;

        amrex::Print() << "Coarse STEP " << step + 1 << " ends."
        << " TIME = " << cur_time << " DT = " << dt_common << std::endl;

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step + 1) % plot_int == 0)
        {
            last_plot_file_step = step + 1;
            plotfilenum++;
            WritePlotFile(plotfilenum);
        }

        if (chk_int > 0 && (step + 1) % chk_int == 0)
        {
            chkfilenum++;
            WriteCheckpointFile(chkfilenum);
        }

        if (monitor_file_int > 0 && (step + 1) % monitor_file_int == 0)
        {
            WriteMonitorFile(cur_time);
        }

        if (cur_time >= stop_time - 1.e-6 * dt_common) break;


        //local cleanup
        flux.clear();
        efield_fc.clear();
        gradne_fc.clear();
        grad_fc.clear();
        expl_src.clear();
        Sborder.clear();
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step)
    {
        plotfilenum++;
        WritePlotFile(plotfilenum);
    }
}
