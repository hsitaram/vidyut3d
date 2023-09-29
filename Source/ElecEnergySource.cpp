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

void Vidyut::compute_elecenergy_source(int lev, const int num_grow, 
                            MultiFab& Sborder, 
                            MultiFab& dsdt,
                            Real time, Real dt)
{
    const auto dx = geom[lev].CellSizeArray();
    auto prob_lo = geom[lev].ProbLoArray();
    auto prob_hi = geom[lev].ProbHiArray();
    ProbParm const* localprobparm = d_prob_parm;

    int ncomp = Sborder.nComp();
    amrex::Real captured_gastemp=gas_temperature;
    amrex::Real captured_gaspres=gas_pressure;

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);

        Array4<Real> sborder_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            //Joule heating
            amrex::Real current_density[AMREX_SPACEDIM];
            amrex::Real Efield[AMREX_SPACEDIM];

            Efield[0]=sborder_arr(i,j,k,EFX_ID);
            Efield[1]=sborder_arr(i,j,k,EFX_ID);
            Efield[2]=sborder_arr(i,j,k,EFX_ID);
             
            //includes sign (negative for electrons/neg-ions)
            amrex::Real mu = plasmachem_transport::mobility(i, j, k, EDN_ID, 
                                  sborder_arr, 
                                  prob_lo, prob_hi, dx, time, 
                                  *localprobparm,captured_gastemp,captured_gaspres);
            
            amrex::Real dcoeff = plasmachem_transport::diffusion_coeff(i, j, k, EDN_ID, 
                                  sborder_arr, 
                                  prob_lo, prob_hi, dx, time, 
                                  *localprobparm,captured_gastemp,captured_gaspres);
            
            amrex::Real nu = plasmachem_transport::collision_freq(i, j, k, EDN_ID, 
                                  sborder_arr, 
                                  prob_lo, prob_hi, dx, time, 
                                  *localprobparm,captured_gastemp,captured_gaspres);


            amrex::Real charge=plasmachem::get_charge(EDN_ID)*ECHARGE;
            amrex::Real ne=sborder_arr(i,j,k,EDN_ID);

            current_density[0]=charge*(mu*ne*Efield[0]-dcoeff*sborder_arr(i,j,k,EDGX_ID));
            current_density[1]=charge*(mu*ne*Efield[1]-dcoeff*sborder_arr(i,j,k,EDGY_ID));
            current_density[2]=charge*(mu*ne*Efield[2]-dcoeff*sborder_arr(i,j,k,EDGZ_ID));

            amrex::Real elec_jheat=current_density[0]*Efield[0]
                        +current_density[1]*Efield[1]+current_density[2]*Efield[2];

            amrex::Real molwt_bg=plasmachem::get_bg_molwt(i, j, k, sborder_arr, *localprobparm);

            //amrex::Real electemp=2.0/3.0*sborder_arr(i,j,k,EEN_ID)/sborder_arr(i,j,k,EDN_ID)/K_B;
            amrex::Real electemp=sborder_arr(i,j,k,ETEMP_ID);

            amrex::Real elec_elastic_coll_term= 3.0/2.0 
                * sborder_arr(i,j,k,EDN_ID) 
                * (electemp-captured_gastemp) * nu;

            //electron energy loss is -ve
            amrex::Real elec_inelastic_coll_term =  plasmachem_reactions::compute_electron_inelastic_heating
                                (i, j, k, 
                                  sborder_arr, 
                                  prob_lo, prob_hi, dx, time, 
                                  *localprobparm,captured_gastemp,
                                  captured_gaspres);

            dsdt_arr(i, j, k) += (elec_jheat - elec_elastic_coll_term + elec_inelastic_coll_term);
               
        });
    }
}
