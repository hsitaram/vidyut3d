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
                            Array<MultiFab,AMREX_SPACEDIM>& efield, 
                            Array<MultiFab,AMREX_SPACEDIM>& gradne, 
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
    
    const int* domlo_arr = geom[lev].Domain().loVect();
    const int* domhi_arr = geom[lev].Domain().hiVect();

    for (MFIter mfi(dsdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx, 1);

        Array4<Real> sborder_arr = Sborder.array(mfi);
        Array4<Real> dsdt_arr = dsdt.array(mfi);
        Array4<Real> phi_arr = phi_new[lev].array(mfi);
        
        GpuArray<Array4<Real>, AMREX_SPACEDIM> 
        ef_arr{AMREX_D_DECL(efield[0].array(mfi), 
                            efield[1].array(mfi), 
                            efield[2].array(mfi))};
        
        GpuArray<Array4<Real>, AMREX_SPACEDIM> 
        gradne_arr{AMREX_D_DECL(gradne[0].array(mfi), 
                                gradne[1].array(mfi), 
                                gradne[2].array(mfi))};

        // update residual
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    
            //Joule heating
            amrex::Real mu,dcoeff,etemp,ne;
            amrex::Real efield_x,efield_y,efield_z,efield_face,gradne_face;
            amrex::Real charge=plasmachem::get_charge(EDN_ID)*ECHARGE;
            amrex::Real current_density[AMREX_SPACEDIM]={0.0};
            amrex::Real efield_fc[AMREX_SPACEDIM]={0.0};
            amrex::Real elec_jheat=0.0;

            //FIXME: This can be done more efficiently sweeping over
            //faces
            for(int idim=0;idim<AMREX_SPACEDIM;idim++)
            {
                //left,right,top,bottom,front,back
                for(int f=0;f<2;f++)
                {
                    IntVect face(i,j,k);
                    face[idim]+=f;

                    if(face[idim]==domlo_arr[idim])
                    {
                       //do the interior face
                       face[idim]+=1;
                    }
                    
                    if(face[idim]==(domhi_arr[idim]+1))
                    {
                       //do the interior face
                        face[idim]-=1;
                    }
                    
                    IntVect lcell(face[0],face[1],face[2]);
                    IntVect rcell(face[0],face[1],face[2]);
                    lcell[idim]-=1;
                    
                    etemp=0.5*(sborder_arr(lcell,ETEMP_ID) 
                               + sborder_arr(rcell,ETEMP_ID));

                    //FIXME:use face centered updated efields here
                    efield_x=0.5*(sborder_arr(lcell,EFX_ID) 
                                  + sborder_arr(rcell,EFX_ID));

                    efield_y=0.5*(sborder_arr(lcell,EFY_ID) 
                                  + sborder_arr(rcell,EFY_ID));

                    efield_z=0.5*(sborder_arr(lcell,EFZ_ID) 
                                  + sborder_arr(rcell,EFZ_ID));

                    ne = 0.5*(sborder_arr(lcell,EDN_ID) 
                              + sborder_arr(rcell,EDN_ID));

                    efield_face=ef_arr[idim](face);
                    gradne_face=gradne_arr[idim](face);

                    mu = plasmachem_transport::mobility(EDN_ID,etemp,
                                                        efield_x,efield_y,efield_z, 
                                                        prob_lo, prob_hi, dx, time, 
                                                        *localprobparm,captured_gastemp,
                                                        captured_gaspres);

                    dcoeff = plasmachem_transport::diffusion_coeff(EDN_ID,etemp,
                                                                   efield_x,efield_y,efield_z, 
                                                                   prob_lo, prob_hi, dx, time, 
                                                                   *localprobparm,captured_gastemp,
                                                                   captured_gaspres);


                    //current_density[idim]+=charge*(mu*ne*efield_face-dcoeff*gradne_face);
                    //efield_fc[idim]+=efield_face;

                    current_density[idim] = charge*(mu*ne*efield_face-dcoeff*gradne_face);
                    efield_fc[idim] = efield_face;

                    elec_jheat += current_density[idim]*efield_face;
                }
            }

            //why only 0.5? 
            //because left and right contribute to JxEx, 
            //top and bottom contribute to JyEy,
            //front and back contribute to JzEz,
            //so we are averaging each component, so it is 1/2 and not 1/6
            //read Deconinck, T, S. Mahadevan, and L. L. Raja. 
            //"Discretization of the Joule heating term for plasma discharge 
            //fluid models in unstructured meshes." 
            //Journal of computational physics 228.12 (2009): 4435-4443.

            elec_jheat*=0.5;
            /*elec_jheat=0.25*(current_density[0]*efield_fc[0]+
                             current_density[1]*efield_fc[1]+
                             current_density[2]*efield_fc[2]);*/
            
            amrex::Real nu = plasmachem_transport::collision_freq(i, j, k, EDN_ID,
                                                                  sborder_arr,
                                                                  prob_lo, prob_hi, dx, time,
                                                                  *localprobparm,captured_gastemp,captured_gaspres);

            amrex::Real molwt_bg=plasmachem::get_bg_molwt(i, j, k, sborder_arr, *localprobparm);

            //amrex::Real electemp=2.0/3.0*sborder_arr(i,j,k,EEN_ID)/sborder_arr(i,j,k,EDN_ID)/K_B;
            amrex::Real electemp=sborder_arr(i,j,k,ETEMP_ID);

            amrex::Real elec_elastic_coll_term= 3.0/2.0 * K_B * ne
            * (electemp-captured_gastemp) * nu * (2.0*ME/molwt_bg);

            //electron energy loss is -ve
            amrex::Real elec_inelastic_coll_term =  plasmachem_reactions::compute_electron_inelastic_heating
            (i, j, k, 
             sborder_arr, 
             prob_lo, prob_hi, dx, time, 
             *localprobparm,captured_gastemp,
             captured_gaspres);

            dsdt_arr(i, j, k) += (elec_jheat - elec_elastic_coll_term + elec_inelastic_coll_term);
            //dsdt_arr(i, j, k) += elec_jheat;
            //
            phi_arr(i,j,k,EJH_ID)=elec_jheat;
            phi_arr(i,j,k,EIH_ID)=elec_inelastic_coll_term;
            phi_arr(i,j,k,EEH_ID)=elec_elastic_coll_term;

        });
    }
}
