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
#include <Transport.H>
#include <Reactions.H>
#include <compute_flux_3d.H>
#include <AMReX_MLABecLaplacian.H>

void echemAMR::compute_dsdt(int lev, const int num_grow, int specid, MultiFab& Sborder, 
                            Array<MultiFab,AMREX_SPACEDIM>& flux, MultiFab& dsdt,
                            Real time, Real dt)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int ncomp = Sborder.nComp();
    int captured_specid = specid;

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);
        FArrayBox reactsource_fab(bx, 1);

        Elixir reactsource_fab_eli = reactsource_fab.elixir();

        Array4<Real> sborder_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);
        Array4<Real> reactsource_arr = reactsource_fab.array();

        GpuArray<Array4<Real>, AMREX_SPACEDIM> flux_arr{
            AMREX_D_DECL(flux[0].array(mfi), 
                         flux[1].array(mfi), flux[2].array(mfi))};

        reactsource_fab.setVal<RunOn::Device>(0.0);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plasmachem_reactions::compute_react_source(i, j, k, captured_specid, sborder_arr, 
                                                       reactsource_arr, prob_lo, prob_hi, dx, time, *localprobparm);
        });

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            dsdt_arr(i, j, k) = (flux_arr[0](i, j, k) - flux_arr[0](i + 1, j, k)) / dx[0] 
            + (flux_arr[1](i, j, k) - flux_arr[1](i, j + 1, k)) / dx[1]
            + (flux_arr[2](i, j, k) - flux_arr[2](i, j, k + 1)) / dx[2] 
            + reactsource_arr(i,j,k);
        });
    }
}

void echemAMR::update_explsrc_at_all_levels(int specid, Vector<MultiFab>& Sborder,
                                            Vector<Array<MultiFab,AMREX_SPACEDIM>>& flux,
                                            Vector<MultiFab>& expl_src, 
                                            amrex::Real cur_time)
{
    for(int lev=0; lev <= finest_level; lev++)
    {
        expl_src[lev].setVal(0.0);
        flux[lev][0].setVal(0.0);
        flux[lev][1].setVal(0.0);
        flux[lev][2].setVal(0.0);
    }

    for(int lev=0; lev <= finest_level; lev++)
    {
        compute_specie_transport_flux(lev, Sborder[lev].nGrow(), Sborder[lev], 
                                      flux[lev], cur_time, specid);
    }

    // =======================================================
    // Average down the fluxes before using them to update phi
    // =======================================================
    for (int lev = finest_level; lev > 0; lev--)
    {
        average_down_faces(amrex::GetArrOfConstPtrs(flux[lev  ]),
                           amrex::GetArrOfPtrs(flux[lev-1]),
                           refRatio(lev-1), Geom(lev-1));
    }

    for(int lev=0;lev<=finest_level;lev++)
    {
        //FIXME: need to avoid this fillpatch
        compute_dsdt(lev, Sborder[lev].nGrow(), specid, Sborder[lev], 
                     flux[lev], expl_src[lev], 
                     cur_time, dt[lev]);
    }
}

void echemAMR::compute_specie_transport_flux(int lev, const int num_grow, MultiFab& Sborder, 
                                             Array<MultiFab,AMREX_SPACEDIM>& flux, 
                                             Real time,int specid)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int ncomp = Sborder.nComp();
    int captured_specid = specid;

    // Get the boundary ids
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();

    GpuArray<int,AMREX_SPACEDIM> domlo={domlo_arr[0], domlo_arr[1], domlo_arr[2]};
    GpuArray<int,AMREX_SPACEDIM> domhi={domhi_arr[0], domhi_arr[1], domhi_arr[2]};
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);
            Box bx_x = convert(bx, IntVect::TheDimensionVector(0));
            Box bx_y = convert(bx, IntVect::TheDimensionVector(1));
            Box bx_z = convert(bx, IntVect::TheDimensionVector(2));

            FArrayBox velx_fab(bx_x, 1);
            FArrayBox vely_fab(bx_y, 1);
            FArrayBox velz_fab(bx_z, 1);

            Elixir velx_fab_eli = velx_fab.elixir();
            Elixir vely_fab_eli = vely_fab.elixir();
            Elixir velz_fab_eli = velz_fab.elixir();

            Array4<Real> sborder_arr = Sborder.array(mfi);
            Array4<Real> velx_arr = velx_fab.array();
            Array4<Real> vely_arr = vely_fab.array();
            Array4<Real> velz_arr = velz_fab.array();

            GpuArray<Array4<Real>, AMREX_SPACEDIM> 
            flux_arr{AMREX_D_DECL(flux[0].array(mfi), 
                                  flux[1].array(mfi), flux[2].array(mfi))};

            velx_fab.setVal<RunOn::Device>(0.0);
            vely_fab.setVal<RunOn::Device>(0.0);
            velz_fab.setVal<RunOn::Device>(0.0);

            amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                plasmachem_transport::compute_vel(i, j, k, 0, captured_specid, sborder_arr, velx_arr, 
                                                  prob_lo, prob_hi, domlo, domhi, dx, time, *localprobparm);
            });

            amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                plasmachem_transport::compute_vel(i, j, k, 1, captured_specid, sborder_arr, vely_arr, 
                                                  prob_lo, prob_hi, domlo, domhi, dx, time, *localprobparm);
            });

            amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                plasmachem_transport::compute_vel(i, j, k, 2, captured_specid, sborder_arr, velz_arr, 
                                                  prob_lo, prob_hi, domlo, domhi, dx, time, *localprobparm);
            });

            amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 0, captured_specid, sborder_arr, 
                             velx_arr, flux_arr[0], dx, *localprobparm); 
            });

            amrex::ParallelFor(bx_y, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 1, captured_specid, sborder_arr, 
                             vely_arr, flux_arr[1], dx, *localprobparm); 
            });

            amrex::ParallelFor(bx_z, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                compute_flux(i, j, k, 2, captured_specid, sborder_arr, 
                             velz_arr, flux_arr[2], dx, *localprobparm);
            });
        }
    }
}

void echemAMR::implicit_solve_species(Real current_time, Real dt, int spec_id, 
                                      Vector<MultiFab>& Sborder, Vector<MultiFab>& dsdt_expl, 
                                      Vector<int>& bc_lo, Vector<int>& bc_hi)
{
    BL_PROFILE("echemAMR::implicit_solve_species(" + std::to_string( spec_id ) + ")");


    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int verbose = 1;
    int captured_spec_id=spec_id;
    int electron_flag=(spec_id==EDN_ID)?1:0;

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
    ProbParm const* localprobparm = d_prob_parm;

    const Real tol_rel = linsolve_reltol;
    const Real tol_abs = linsolve_abstol;

    // set A and B, A=1/dt, B=1
    Real ascalar = 1.0/dt;
    Real bscalar = 1.0;

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_lo 
    = {LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_linsolve_hi 
    = {LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin}; 

    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (bc_lo[idim] == PERBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (bc_lo[idim] == DIRCBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_lo[idim] == HNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Neumann;
        }
        if (bc_lo[idim] == IHNEUBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_lo[idim] == ROBINBC)
        {
            bc_linsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }

        //higher side bcs
        if (bc_hi[idim] == PERBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Periodic;
        }
        if (bc_hi[idim] == DIRCBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_hi[idim] == HNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Neumann;
        }
        if (bc_hi[idim] == IHNEUBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_hi[idim] == ROBINBC)
        {
            bc_linsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
    }

    Vector<MultiFab> specdata(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<Array<MultiFab, AMREX_SPACEDIM>> gradsoln(finest_level+1);
    Vector<MultiFab> bcoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        specdata[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        bcoeff[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, 1);

        if(electron_flag)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& faceba = amrex::convert(grids[ilev], 
                                                        IntVect::TheDimensionVector(idim));
                gradsoln[ilev][idim].define(faceba, dmap[ilev], 1, 0);
            }
        }
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    linsolve_ptr.reset(new MLABecLaplacian(Geom(0,finest_level), 
                                           boxArray(0,finest_level), 
                                           DistributionMap(0,finest_level), info));
    MLMG mlmg(*linsolve_ptr);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(verbose);
    linsolve_ptr->setDomainBC(bc_linsolve_lo, bc_linsolve_hi);
    linsolve_ptr->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        // Copy args (FabArray<FAB>& dst, FabArray<FAB> const& src, 
        // int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
        specdata[ilev].setVal(0.0);
        amrex::Copy(specdata[ilev], Sborder[ilev], captured_spec_id, 0, 1, num_grow);

        acoeff[ilev].setVal(1.0);
        bcoeff[ilev].setVal(1.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        rhs[ilev].setVal(0.0);
        MultiFab::LinComb(rhs[ilev], 1.0/dt, specdata[ilev], 0, 1.0, 
                          dsdt_expl[ilev], 0, 0, 1, 0);

        solution[ilev].setVal(0.0);
        amrex::MultiFab::Copy(solution[ilev], specdata[ilev], 0, 0, 1, 0);
        int ncomp = Sborder[ilev].nComp();

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

            Array4<Real> phi_arr = Sborder[ilev].array(mfi);
            Array4<Real> acoeff_arr = acoeff[ilev].array(mfi);
            Array4<Real> bcoeff_arr = bcoeff[ilev].array(mfi);

            FArrayBox dcoeff_fab(gbx, 1);
            Elixir dcoeff_fab_eli = dcoeff_fab.elixir();
            Array4<Real> dcoeff_arr = dcoeff_fab.array();

            dcoeff_fab.setVal<RunOn::Device>(0.0);

            amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                //FIXME: use component wise call
                plasmachem_transport::compute_dcoeff(i, j, k, captured_spec_id, phi_arr, 
                                                     dcoeff_arr, prob_lo, 
                                                     prob_hi, dx, time, *localprobparm);

                bcoeff_arr(i,j,k)=dcoeff_arr(i,j,k);
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

            Array4<Real> bc_arr = specdata[ilev].array(mfi);
            Array4<Real> phi_arr = Sborder[ilev].array(mfi);
            Real time = current_time; // for GPU capture

            Array4<Real> robin_a_arr = robin_a[ilev].array(mfi);
            Array4<Real> robin_b_arr = robin_b[ilev].array(mfi);
            Array4<Real> robin_f_arr = robin_f[ilev].array(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (!geom[ilev].isPeriodic(idim))
                {
                    //note: bdryLo/bdryHi grabs the face indices from bx that are the boundary
                    //since they are face indices, the bdry normal index is 0/n+1, n is number of cells
                    //so the ghost cell index at left side is i-1 while it is i on the right
                    if (bx.smallEnd(idim) == domain.smallEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryLo(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            plasmachem_transport::species_bc(i, j, k, idim, -1, 
                                                             captured_spec_id, phi_arr, bc_arr, robin_a_arr,
                                                             robin_b_arr, robin_f_arr, 
                                                             prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                    if (bx.bigEnd(idim) == domain.bigEnd(idim))
                    {
                        amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            plasmachem_transport::species_bc(i, j, k, idim, +1, 
                                                             captured_spec_id, phi_arr, bc_arr, robin_a_arr, 
                                                             robin_b_arr, robin_f_arr,
                                                             prob_lo, prob_hi, dx, time, *localprobparm);
                        });
                    }
                }
            }
        }

        linsolve_ptr->setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_ptr->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells
        if(mixedbc)
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]), &(robin_a[ilev]), &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_ptr->setLevelBC(ilev, &(specdata[ilev]));
        }
    }

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    //bound species density
    if(bound_specden)
    { 
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            amrex::Real minspecden=min_species_density; 
            for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> soln_arr = solution[ilev].array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    if(soln_arr(i,j,k) < minspecden)
                    {
                        soln_arr(i,j,k)=minspecden;
                    }
                });
            }
        }
    }
    if(electron_flag)
    {
        mlmg.getGradSolution(GetVecOfArrOfPtrs(gradsoln));
    }

    // copy solution back to phi_new
    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, spec_id, 1, 0);
    }
    Print()<<"Solved species:"<<allvarnames[spec_id]<<"\n";

    if(electron_flag)
    {
        for (int ilev = 0; ilev <= finest_level; ilev++)
        {
            const Array<const MultiFab*, AMREX_SPACEDIM> allgrad = {&gradsoln[ilev][0], 
                &gradsoln[ilev][1], &gradsoln[ilev][2]};
            average_face_to_cellcenter(phi_new[ilev], EDGX_ID, allgrad);
        }
    }

    //clean-up
    specdata.clear();
    acoeff.clear();
    bcoeff.clear();
    gradsoln.clear();
    solution.clear();
    rhs.clear();

    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}

