#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLTensorOp.H>

#include <Vidyut.H>
#include <UserSources.H>
#include <Chemistry.H>
#include <UserFunctions.H>
#include <PlasmaChem.H>
#include <ProbParm.H>
#include <AMReX_MLABecLaplacian.H>

void Vidyut::solve_potential(Real current_time, Vector<MultiFab>& Sborder,
                             amrex::Vector<int>& bc_lo,amrex::Vector<int>& bc_hi,
                             Vector<Array<MultiFab,AMREX_SPACEDIM>>& efield_fc)
{
    BL_PROFILE("Vidyut::solve_potential()");

    // FIXME: add these as inputs
    int max_coarsening_level = linsolve_max_coarsening_level;
    int max_iter=linsolve_maxiter;
    Real ascalar = 1.0;
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
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;

    // default to inhomogNeumann since it is defaulted to flux = 0.0 anyways
    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_lo 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    std::array<LinOpBCType, AMREX_SPACEDIM> bc_potsolve_hi 
    = {AMREX_D_DECL(LinOpBCType::Robin, LinOpBCType::Robin, LinOpBCType::Robin)}; 

    int mixedbc=0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
        //lower side bcs
        if (bc_lo[idim] == PERBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Periodic;
        }
        if (bc_lo[idim] == DIRCBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_lo[idim] == HNEUBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Neumann;
        }
        if (bc_lo[idim] == IHNEUBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_lo[idim] == ROBINBC)
        {
            bc_potsolve_lo[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }

        //higher side bcs
        if (bc_hi[idim] == PERBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Periodic;
        }
        if (bc_hi[idim] == DIRCBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Dirichlet;
        }
        if (bc_hi[idim] == HNEUBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Neumann;
        }
        if (bc_hi[idim] == IHNEUBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::inhomogNeumann;
        }
        if (bc_hi[idim] == ROBINBC)
        {
            bc_potsolve_hi[idim] = LinOpBCType::Robin;
            mixedbc=1;
        }
    }

    Vector<MultiFab> potential(finest_level+1);
    Vector<MultiFab> acoeff(finest_level+1);
    Vector<MultiFab> solution(finest_level+1);
    Vector<MultiFab> rhs(finest_level+1);

    Vector<MultiFab> robin_a(finest_level+1);
    Vector<MultiFab> robin_b(finest_level+1);
    Vector<MultiFab> robin_f(finest_level+1);

    const int num_grow = 1;

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        potential[ilev].define(grids[ilev], dmap[ilev], 1, num_grow);
        acoeff[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        robin_a[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_b[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        robin_f[ilev].define(grids[ilev], dmap[ilev], 1, 1);
    }

    LPInfo info;
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMaxCoarseningLevel(max_coarsening_level);
    linsolve_ptr.reset(new MLABecLaplacian(Geom(0,finest_level), 
                                           boxArray(0,finest_level), 
                                           DistributionMap(0,finest_level), info));

    linsolve_ptr->setMaxOrder(2);
    linsolve_ptr->setDomainBC(bc_potsolve_lo, bc_potsolve_hi);
    linsolve_ptr->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        potential[ilev].setVal(0.0);

        // Copy (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, 
        // int dstcomp, int numcomp, const IntVect& nghost)
        amrex::Copy(potential[ilev], Sborder[ilev], POT_ID, 0, 1, num_grow);

        solution[ilev].setVal(0.0);
        // FIXME: for some reason copying in current soln breaks the solver...
        // amrex::MultiFab::Copy(solution[ilev], potential[ilev], 0, 0, 1, 0);
        rhs[ilev].setVal(0.0);
        acoeff[ilev].setVal(0.0);

        //default to homogenous Neumann
        robin_a[ilev].setVal(0.0);
        robin_b[ilev].setVal(1.0);
        robin_f[ilev].setVal(0.0);

        // Get the boundary ids
        const int* domlo_arr = geom[ilev].Domain().loVect();
        const int* domhi_arr = geom[ilev].Domain().hiVect();

        GpuArray<int,AMREX_SPACEDIM> domlo={AMREX_D_DECL(domlo_arr[0], domlo_arr[1], domlo_arr[2])};
        GpuArray<int,AMREX_SPACEDIM> domhi={AMREX_D_DECL(domhi_arr[0], domhi_arr[1], domhi_arr[2])};


        // fill cell centered diffusion coefficients and rhs
        if(do_spacechrg){
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
                Array4<Real> rhs_arr = rhs[ilev].array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    //include electrons
                    rhs_arr(i,j,k)=0.0;
                    for(int sp=0;sp<NUM_SPECIES;sp++)
                    {
                        if(amrex::Math::abs(plasmachem::get_charge(sp)) > 0)
                        {
                            rhs_arr(i,j,k)+=plasmachem::get_charge(sp)*phi_arr(i,j,k,sp);
                        }
                    }
                    //why minus sign?
                    //remember del.E = rho/eps0
                    //but E= -grad phi
                    //del^2 phi = -rho/eps0
                    rhs_arr(i,j,k)*=(-ECHARGE/EPS0);

                    user_sources::add_user_potential_sources(i, j, k, phi_arr, 
                                                             rhs_arr, prob_lo, prob_hi, 
                                                             dx, time, *localprobparm);
                });
            }
        }

        // average cell coefficients to faces, this includes boundary faces
        Array<MultiFab, AMREX_SPACEDIM> face_bcoeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(acoeff[ilev].boxArray(), 
                                                IntVect::TheDimensionVector(idim));
            face_bcoeff[idim].define(ba, acoeff[ilev].DistributionMap(), 1, 0);
            face_bcoeff[idim].setVal(-1.0); //since bscalar was set to 1
        }
        // set boundary conditions
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto dx = geom[ilev].CellSizeArray();
            auto prob_lo = geom[ilev].ProbLoArray();
            auto prob_hi = geom[ilev].ProbHiArray();
            const Box& domain = geom[ilev].Domain();

            Array4<Real> phi_arr = Sborder[ilev].array(mfi);
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
                        int domend = -1;
                        amrex::Real app_voltage = get_applied_potential(time, domend);
                        if(user_defined_potential == 1){
                            user_transport::potential_bc(i, j, k, idim, -1, 
                                                         phi_arr, bc_arr, robin_a_arr, 
                                                         robin_b_arr, robin_f_arr, 
                                                         prob_lo, prob_hi, dx, time, 
                                                         *localprobparm,captured_gastemp,
                                                         captured_gaspres, app_voltage);
                        } else {
                            plasmachem_transport::potential_bc(i, j, k, idim, -1, 
                                                               phi_arr, bc_arr, robin_a_arr, 
                                                               robin_b_arr, robin_f_arr, 
                                                               prob_lo, prob_hi, dx, time, 
                                                               *localprobparm,captured_gastemp,
                                                               captured_gaspres, app_voltage);
                        }
                    });
                }
                if (bx.bigEnd(idim) == domain.bigEnd(idim))
                {
                    amrex::ParallelFor(amrex::bdryHi(bx, idim), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        int domend = 1;
                        amrex::Real app_voltage = get_applied_potential(time, domend);
                        if(user_defined_potential == 1){
                            user_transport::potential_bc(i, j, k, idim, +1, 
                                                         phi_arr, bc_arr, robin_a_arr, 
                                                         robin_b_arr, robin_f_arr, 
                                                         prob_lo, prob_hi, dx, time,
                                                         *localprobparm,captured_gastemp,
                                                         captured_gaspres, app_voltage);
                        } else {
                            plasmachem_transport::potential_bc(i, j, k, idim, +1, 
                                                               phi_arr, bc_arr, robin_a_arr, 
                                                               robin_b_arr, robin_f_arr, 
                                                               prob_lo, prob_hi, dx, time,
                                                               *localprobparm,captured_gastemp,
                                                               captured_gaspres, app_voltage);
                        }
                    });
                }
            }
        }

        linsolve_ptr->setACoeffs(ilev, acoeff[ilev]);

        // set b with diffusivities
        linsolve_ptr->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoeff));

        // bc's are stored in the ghost cells of potential
        if(mixedbc)
        {
            linsolve_ptr->setLevelBC(ilev, &potential[ilev], &(robin_a[ilev]), 
                                     &(robin_b[ilev]), &(robin_f[ilev]));
        }
        else
        {
            linsolve_ptr->setLevelBC(ilev, &potential[ilev]);
        }

    }

    MLMG mlmg(*linsolve_ptr);
    mlmg.setMaxIter(linsolve_maxiter);
    mlmg.setVerbose(verbose);
    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    mlmg.getGradSolution(GetVecOfArrOfPtrs(efield_fc));

    amrex::Print()<<"Solved Potential\n";

    for (int ilev = 0; ilev <= finest_level; ilev++)
    {
        amrex::MultiFab::Copy(phi_new[ilev], solution[ilev], 0, POT_ID, 1, 0);
        efield_fc[ilev][0].mult(-1.0,0,1);
#if AMREX_SPACEDIM > 1
        efield_fc[ilev][1].mult(-1.0,0,1);
#if AMREX_SPACEDIM == 3
        efield_fc[ilev][2].mult(-1.0,0,1);
#endif
#endif   

        const Array<const MultiFab*, AMREX_SPACEDIM> allgrad = {AMREX_D_DECL(&efield_fc[ilev][0], 
            &efield_fc[ilev][1], &efield_fc[ilev][2])};
        average_face_to_cellcenter(phi_new[ilev], EFX_ID, allgrad);

        // Calculate the reduced electric field
        for (MFIter mfi(phi_new[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, 1);

            Array4<Real> s_arr = phi_new[ilev].array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                RealVect Evect{AMREX_D_DECL(s_arr(i,j,k,EFX_ID),s_arr(i,j,k,EFY_ID),s_arr(i,j,k,EFZ_ID))};
                Real Esum = 0.0;
                for(int dim=0; dim<AMREX_SPACEDIM; dim++) Esum += Evect[dim]*Evect[dim];
                s_arr(i,j,k,REF_ID) = (pow(Esum, 0.5) / gas_num_dens) / 1.0e-21;
            });
        }
    }

    //clean-up
    potential.clear();
    acoeff.clear();
    solution.clear();
    rhs.clear();

    robin_a.clear();
    robin_b.clear();
    robin_f.clear();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Vidyut::get_applied_potential(Real current_time, int domain_end)
{
    ProbParm const* localprobparm = d_prob_parm;
    Real voltage;
    Real voltage_amp = (domain_end == -1) ? voltage_amp_1:voltage_amp_2;
    if(voltage_profile == 1) {  // Sinusoidal pulse shape
        voltage = sin(2.0*PI*voltage_freq*current_time)*voltage_amp;
    } else if (voltage_profile == 2) {    // Single triangular pulse
        if(current_time <= voltage_center) {
            voltage = (voltage_center - current_time < voltage_dur/2.0) ? (1.0 - (voltage_center - current_time)/(voltage_dur/2.0))*voltage_amp:0.0;
        } else {
            voltage = (current_time - voltage_center < voltage_dur/2.0) ? (1.0 - (current_time - voltage_center)/(voltage_dur/2.0))*voltage_amp:0.0;
        }
    } else {
        voltage = (domain_end==-1) ? voltage_amp_1:voltage_amp_2;
    }

    return voltage;
}
