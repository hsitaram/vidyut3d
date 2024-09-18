#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <global_defines.H>

using namespace amrex;

void initlevset(MultiFab &levset,Geometry geom)
{
    Real int_loc=0.5;
    Real rad1=0.2;
    Real rad2=0.4;
    int transport_dir=0;

    ParmParse pp;
    pp.query("rad1",rad1);
    pp.query("rad2",rad2);

    for (MFIter mfi(levset); mfi.isValid(); ++mfi)
    {
        const auto dx = geom.CellSizeArray();
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx,1);
        auto prblo = geom.ProbLoArray();
        auto prbhi = geom.ProbHiArray();

        Array4<Real> ls_arr = levset.array(mfi);

        amrex::ParallelFor(gbx,[=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real midx=0.5*(prblo[0]+prbhi[0]);
            Real midy=0.5*(prblo[1]+prbhi[1]);
            Real midz=0.5*(prblo[2]+prbhi[2]);

            Real x[AMREX_SPACEDIM];
            Real spltcell[AMREX_SPACEDIM];

            //split cell approach
            /*spltcell[0]=8;
              spltcell[1]=8;
              spltcell[2]=8;

              Real vfrac=0.0;

              for (int kk=0;kk<spltcell[2];kk++)
              {
              for(int jj=0;jj<spltcell[1];jj++)
              {
              for(int ii=0;ii<spltcell[0];ii++)
              {
              x[0]= prblo[0]+(i+(ii+0.5)/spltcell[0])*dx[0];
              x[1]= prblo[1]+(j+(jj+0.5)/spltcell[1])*dx[1];
              x[2]= prblo[2]+(k+(kk+0.5)/spltcell[2])*dx[2];


              Real dist=std::sqrt((x[0]-midx)*(x[0]-midx)+
              (x[1]-midy)*(x[1]-midy)+
              (x[2]-midz)*(x[2]-midz));


              if(dist>rad1 && dist<rad2)
              {
              vfrac+=1.0;
              }
              }
              }
              }*/
            //ls_arr(i,j,k) = vfrac/Real(spltcell[0]*spltcell[1]*spltcell[2]);

            //corner approach
            Real vfrac=0.0;

            for (int kk=0;kk<2;kk++)
            {
                for(int jj=0;jj<2;jj++)
                {
                    for(int ii=0;ii<2;ii++)
                    {
                        x[0]= prblo[0]+(i+ii)*dx[0];
                        x[1]= prblo[1]+(j+jj)*dx[1];
                        x[2]= prblo[2]+(k+kk)*dx[2];


                        Real dist=std::sqrt((x[0]-midx)*(x[0]-midx)+
                                            (x[1]-midy)*(x[1]-midy)+
                                            (x[2]-midz)*(x[2]-midz));


                        if(dist>rad1 && dist<rad2)
                        {
                            vfrac+=1.0;
                        }
                    }
                }
            }

            ls_arr(i,j,k) = 0.125*vfrac;
        });
    }
}

void setbc(MultiFab &phi_bc,MultiFab &phi, Geometry geom,
           amrex::Vector<int> bc_lo,amrex::Vector<int> bc_hi)
{
    ParmParse pp;
    //pp.query("electrode_dcoeff",dcoeffleft);

    amrex::MultiFab::Copy(phi_bc,phi,0,0,1,0);

    const int* domloarr = geom.Domain().loVect();
    const int* domhiarr = geom.Domain().hiVect();

    GpuArray<int,AMREX_SPACEDIM> domlo={domloarr[0],domloarr[1],domloarr[2]};
    GpuArray<int,AMREX_SPACEDIM> domhi={domhiarr[0],domhiarr[1],domhiarr[2]};

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const auto dx = geom.CellSizeArray();
        const Box& bx = mfi.tilebox();
        const Box& gbx = amrex::grow(bx,1);
        auto prblo = geom.ProbLoArray();
        auto prbhi = geom.ProbHiArray();
        const Box& domain = geom.Domain();

        Array4<Real> phibc_arr=phi_bc.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k){
            IntVect cellid(i,j,k);
            IntVect gcellid(i,j,k);

            if(cellid[XDIR]==domlo[XDIR])
            {
                gcellid[XDIR] -=1;
                //inhom neumann
                //phibc_arr(gcellid) = -current_input/dcoeffleft;
            }
            if(cellid[XDIR]==domhi[XDIR])
            {
                gcellid[XDIR] += 1;
                phibc_arr(gcellid) = 0.0; //dirichlet
            }

            if(cellid[YDIR]==domlo[YDIR])
            {
                gcellid[YDIR] -=1;
                if(bc_lo[YDIR]==DIRICHLET_ID)
                {
                    phibc_arr(gcellid) = 0.0;
                }
            }
            if(cellid[YDIR]==domhi[YDIR])
            {
                gcellid[YDIR] += 1;
                if(bc_hi[YDIR]==DIRICHLET_ID)
                {
                    phibc_arr(gcellid) = 0.0;
                }
            }

            if(cellid[ZDIR]==domlo[ZDIR])
            {
                gcellid[ZDIR] -=1;
                if(bc_lo[ZDIR]==DIRICHLET_ID)
                {
                    phibc_arr(gcellid) = 0.0;
                }
            }
            if(cellid[ZDIR]==domhi[ZDIR])
            {
                gcellid[ZDIR] += 1;
                if(bc_hi[ZDIR]==DIRICHLET_ID)
                {
                    phibc_arr(gcellid) = 0.0;
                }
            }
        });
    }
}
