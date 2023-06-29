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

void echemAMR::solve_potential(Real current_time)
{
    BL_PROFILE("echemAMR::solve_potential()");

    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int max_iter=200;
    Real ascalar = 0.0;
    Real bscalar = 1.0;
    ProbParm const* localprobparm = d_prob_parm;

    //==================================================
    // amrex solves
    // read small a as alpha, b as beta

    //(A a - B del.(b del)) phi = f
    //
    // A and B are scalar constants
    // a and b are scalar fields
    // f is rhs
    // in this case: A=0,a=0,B=1,b=conductivity
    // note also the negative sign
    //====================================================

    // FIXME: had to adjust for constant coefficent,
    // this could be due to missing terms in the 
    // intercalation reaction or sign mistakes...
    const Real tol_rel = linsolve_reltol;
    const Real tol_abs = linsolve_abstol;

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_lo 
    = {LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin};

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_hi 
    = {LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin};

    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        if (bc_lo[idim] == BCType::int_dir)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (bc_hi[idim] == BCType::int_dir)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Periodic;
        }
    }

    Vector<MultiFab> potential;
    Vector<MultiFab> acoeff;
    Vector<MultiFab> bcoeff;
    Vector<Array<MultiFab*, AMREX_SPACEDIM>> gradsoln;
    Vector<MultiFab> solution;
    Vector<MultiFab> rhs;

    Vector<MultiFab> robin_a;
    Vector<MultiFab> robin_b;
    Vector<MultiFab> robin_f;

    acoeff.resize(finest_level + 1);
    bcoeff.resize(finest_level + 1);
    gradsoln.resize(finest_level + 1);
    potential.resize(finest_level + 1);
    solution.resize(finest_level + 1);
    rhs.resize(finest_level + 1);

    robin_a.resize(finest_level+1);
    robin_b.resize(finest_level+1);
    robin_f.resize(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        potential[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        bcoeff[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& faceba = amrex::convert(grids[ilev], 
                                                    IntVect::TheDimensionVector(idim));
            gradsoln[ilev][idim] = new MultiFab(faceba, dmap[ilev], 1, 0);
        }

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, 1);
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    MLABecLaplacian mlabec(geom, grids, dmap, info);
    MLMG mlmg(mlabec);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(verbose);
    mlabec.setDomainBC(bc_potsolve_lo, bc_potsolve_hi);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        MultiFab Sborder(grids[ilev], dmap[ilev], phi_new[ilev].nComp(), num_grow);
        FillPatch(ilev, current_time, Sborder, 0, Sborder.nComp());
        potential[ilev].setVal(0.0);

        // Copy (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, 
        // int dstcomp, int numcomp, const IntVect& nghost)
        amrex::Copy(potential[ilev], Sborder, POT_ID, 0, 1, num_grow);


        //dcoeff is 1
        bcoeff[ilev].setVal(1.0);

        solution[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        // Get the boundary ids
        const int* domlo_arr = geom[ilev].Domain().loVect();
        const int* domhi_arr = geom[ilev].Domain().hiVect();

        GpuArray<int,AMREX_SPACEDIM> domlo={domlo_arr[0], domlo_arr[1], domlo_arr[2]};
        GpuArray<int,AMREX_SPACEDIM> domhi={domhi_arr[0], domhi_arr[1], domhi_arr[2]};


        // fill cell centered diffusion coefficients and rhs
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Real time = current_time; // for GPU capture

            Array4<Real> phi_arr = Sborder.array(mfi);
            Array4<Real> bcoeff_arr = bcoeff[ilev].array(mfi);
            Array4<Real> rhs_arr = rhs[ilev].array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                plasmachem_reactions::compute_potential_source(i, j, k, phi_arr, 
                             rhs_arr, prob_lo, prob_hi, dx, time, *localprobparm);
            });
        }

        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoeff[ilev].boxArray(), IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, bcoeff[ilev].DistributionMap(), 1, 0);
        }
        // true argument for harmonic averaging
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoeff), bcoeff[ilev], geom[ilev], true);

        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> phi_arr = Sborder.array(mfi);
            Array4<Real> bc_arr = potential[ilev].array(mfi);

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            Real time = current_time; // for GPU capture

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                //so the ghost cell index at left side is i-1 while it is i on the right
                if (bx.smallEnd(idim) == domain.smallEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        plasmachem_transport::potential_bc(i, j, k, idim, -1, phi_arr, robin_a_arr, 
                                                           robin_b_arr, robin_f_arr, prob_lo, prob_hi, dx, time, *localprobparm);
                    });
                }
                if (bx.bigEnd(idim) == domain.bigEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        plasmachem_transport::potential_bc(i, j, k, idim, +1, phi_arr, robin_a_arr, 
                                                           robin_b_arr, robin_f_arr, prob_lo, prob_hi, dx, time, *localprobparm);
                    });
                }
            }
        }

        // bc's are stored in the ghost cells of potential
        mlabec.setLevelBC(ilev, &potential[ilev], &(robin_a[ilev]), &(robin_b[ilev]), &(robin_f[ilev]));

        acoeff[ilev].setVal(1.0); //will be scaled by ascalar
        mlabec.setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));
    }
    mlabec.setScalars(ascalar, bscalar);

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    mlmg.getGradSolution(gradsoln);

    // copy solution back to phi_new
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, POT_ID, 1, 0);
        const Array<const MultiFab*, AMREX_SPACEDIM> allgrad = {gradsoln[ilev][0], gradsoln[ilev][1], gradsoln[ilev][2]};
        average_face_to_cellcenter(phi_new[ilev], EFX_ID, allgrad);
        phi_new[ilev].mult(-1.0, EFX_ID, 3);
    }
}
