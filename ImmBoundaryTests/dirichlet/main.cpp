#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>
#include <global_defines.H>
#include <Prob.H>
#include <gmres.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
amrex::Real euclid_dist(amrex::Real x[AMREX_SPACEDIM], amrex::Real xi[AMREX_SPACEDIM])
{
    amrex::Real r;
    r=(x[0]-xi[0])*(x[0]-xi[0]);
#if AMREX_SPACEDIM > 1
    r+=(x[1]-xi[1])*(x[1]-xi[1]);
#if AMREX_SPACEDIM == 3
    r+=(x[2]-xi[2])*(x[2]-xi[2]);
#endif
#endif
   return(std::sqrt(r));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
amrex::Real rbf(amrex::Real r,
                amrex::Real eps)
{
   return(std::exp(-eps*eps*r*r));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE 
amrex::Real rbfder(amrex::Real r, amrex::Real eps, amrex::Real x, amrex::Real x0)
{
    amrex::Real derval=-std::exp(-eps*eps*r*r)*(eps*eps)*2.0*(x-x0);
    return(derval);
}

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        std::vector<int> n_cell;
        std::vector<Real> prob_lo;
        std::vector<Real> prob_hi;
        amrex::Real dircval_rad1=1.0;
        amrex::Real dircval_rad2=0.0;
        amrex::Real sourceterm=0.0;
        
        int max_grid_size = 32;
        int use_hypre=0;
        // typically a huge number so MG coarsens as much as possible
        int max_coarsening_level = 20;    
        int max_iter=200;
        Real penalty=1e4;
        Real rad1=0.2;
        Real rad2=0.4;
        int niter=3;

        Vector<int> is_periodic(AMREX_SPACEDIM,0);
        amrex::Vector<int> bc_lo{0,0,0};
        amrex::Vector<int> bc_hi{0,0,0};

        ParmParse pp;
        pp.getarr("n_cell", n_cell);
        pp.getarr("prob_lo", prob_lo);
        pp.getarr("prob_hi", prob_hi);
        pp.getarr("is_periodic",is_periodic);
        pp.query("rad1",rad1);
        pp.query("rad2",rad2);
        pp.query("dircval_rad1",dircval_rad1);
        pp.query("dircval_rad2",dircval_rad2);
        pp.query("sourceterm",sourceterm);
        pp.query("niter",niter);

        pp.query("max_grid_size", max_grid_size);
        pp.query("max_coarsening_level",max_coarsening_level);
        pp.query("max_iter",max_iter);
        pp.query("penalty",penalty);

        pp.queryarr("lo_bc", bc_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("hi_bc", bc_hi, 0, AMREX_SPACEDIM);

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        RealBox rb({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])}, 
                   {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

        Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                   IntVect{AMREX_D_DECL(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1)});

        geom.define(domain, &rb, CoordSys::cartesian, is_periodic.data());

        grids.define(domain); 
        grids.maxSize(max_grid_size); 
        dmap.define(grids); 
        int required_coarsening_level = 0; 

        MultiFab phi(grids, dmap, 1, 1);
        MultiFab soln(grids, dmap, 1, 1);
        MultiFab rhs(grids, dmap, 1, 1);
        MultiFab phi_bc(grids, dmap, 1, 1);
        MultiFab levset(grids,dmap, 1, 1);
        MultiFab acoef(grids, dmap, 1, 1);
        MultiFab cell_bcoef(grids, dmap, 1, 1);
        MultiFab err(grids, dmap, 1,  0);


        initlevset(levset,geom);
        levset.FillBoundary(geom.periodicity());

        LPInfo info;
        info.setMaxCoarseningLevel(max_coarsening_level);
        MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);
        MLMG mlmg(mlabec);
        mlmg.setMaxIter(max_iter);

#ifdef AMREX_USE_HYPRE
        if(use_hypre)
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            Hypre::Interface hypre_interface = Hypre::Interface::ij;
            amrex::Print() << "using hypre" << std::endl;
            mlmg.setHypreOptionsNamespace("hypre");
            mlmg.setHypreInterface(hypre_interface);
        }
#endif


        // relative and absolute tolerances for linear solve
        const Real tol_rel = 1.e-10;
        const Real tol_abs = 1e-10;
        int verbose = 1;
        mlmg.setVerbose(verbose);

        // define array of LinOpBCType for domain boundary conditions
        std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_bc_lo;
        std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_bc_hi;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
        {
            if(bc_lo[idim] == PERIODIC_ID)
            {
                mlmg_bc_lo[idim]=LinOpBCType::Periodic;
            }
            else if(bc_lo[idim] == DIRICHLET_ID)
            {
                mlmg_bc_lo[idim]=LinOpBCType::Dirichlet;
            }
            else if(bc_lo[idim] == H_NEUMANN_ID)
            {
                mlmg_bc_lo[idim]=LinOpBCType::Neumann;
            }
            else
            {
                mlmg_bc_lo[idim]=LinOpBCType::inhomogNeumann;
            }

            if(bc_hi[idim] == PERIODIC_ID)
            {
                mlmg_bc_hi[idim]=LinOpBCType::Periodic;
            }
            else if(bc_hi[idim] == DIRICHLET_ID)
            {
                mlmg_bc_hi[idim]=LinOpBCType::Dirichlet;
            }
            else if(bc_hi[idim] == H_NEUMANN_ID)
            {
                mlmg_bc_hi[idim]=LinOpBCType::Neumann;
            }
            else
            {
                mlmg_bc_hi[idim]=LinOpBCType::inhomogNeumann;
            }
            
        }
        
        Real rel_errnorm;
        Real errnorm_1st=1.0;
        phi.setVal(0.0);
       
        { 
            MultiFab plotfile_mf(grids, dmap, 2, 0);
            MultiFab::Copy(plotfile_mf, levset,0,0,1,0);
            MultiFab::Copy(plotfile_mf, phi,0,1,1,0);

            WriteSingleLevelPlotfile("plt0", plotfile_mf, {"levset","potential"}, geom, 0.0, 0);
        }
        for(int iter=0;iter<niter;iter++)
        {
            rhs.setVal(0.0);
            soln.setVal(0.0);
            phi_bc.setVal(0.0);
            cell_bcoef.setVal(1.0);
            acoef.setVal(1.0);

            phi.FillBoundary(geom.periodicity());
            amrex::MultiFab::Copy(err,phi, 0, 0, 1 ,0);
            //setbc(phi_bc,phi,geom,bc_lo,bc_hi);


            // Boundary of the whole domain. This functions must be called,
            // and must be called before other bc functions.
            mlabec.setDomainBC(mlmg_bc_lo,mlmg_bc_hi);
            mlabec.setLevelBC(0, &phi_bc);

            // operator looks like (ACoef - div BCoef grad) phi = rhs
            // scaling factors; these multiply ACoef and BCoef
            Real ascalar = 1.0;
            Real bscalar = 1.0;
            mlabec.setScalars(ascalar, bscalar);


            // set BCoef to 1.0 (and array of face-centered coefficients)
            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
            {
                face_bcoef[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0);
                face_bcoef[idim].setVal(1.0);
            }

            for (MFIter mfi(phi); mfi.isValid(); ++mfi)
            {
                const auto dx = geom.CellSizeArray();
                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real> phi_arr = phi.array(mfi);
                Array4<Real> ls_arr = levset.array(mfi);
                auto prblo = geom.ProbLoArray();
                auto prbhi = geom.ProbHiArray();

                amrex::ParallelFor(gbx,[=] AMREX_GPU_DEVICE (int i, int j, int k) {

                    if(ls_arr(i,j,k)>0.0 && ls_arr(i,j,k)<1.0)
                    {
                        //amrex::Print()<<"ls,i,j,k:"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<ls_arr(i,j,k)<<"\n";
                        Real midx[AMREX_SPACEDIM];
                        midx[0]=0.5*(prblo[0]+prbhi[0]);
                        midx[1]=0.5*(prblo[1]+prbhi[1]);
                        midx[2]=0.5*(prblo[2]+prbhi[2]);

                        Real ccx[AMREX_SPACEDIM];

                        ccx[0]=prblo[0]+(i+0.5)*dx[0];
                        ccx[1]=prblo[1]+(j+0.5)*dx[1];
                        ccx[2]=prblo[2]+(k+0.5)*dx[2];

                        Real dist=euclid_dist(ccx,midx);

                        Real rad1dist=rad1-dist;
                        Real rad2dist=dist-rad2;

                        Real tol=1e-12;
                        Real mindist=(amrex::Math::abs(rad1dist)<amrex::Math::abs(rad2dist))?rad1dist:rad2dist;
                        bool inner_radflag=(mindist==rad1dist)?true:false;

                        Real mat[NMAX][NMAX]={{0.0}};
                        Real bvec[NMAX]={0.0};
                        Real soln0[NMAX]={0.0};
                        Real wts[NMAX]={0.0};
                        Real eps=1.0/(std::min(dx[0],std::min(dx[1],dx[2])));
                        int stenc_i[NMAX]={0};
                        int stenc_j[NMAX]={0};
                        int stenc_k[NMAX]={0};

                        //normal along r here
                        Real nx=(ccx[0]-midx[0])/dist;
                        Real ny=(ccx[1]-midx[1])/dist;
                        Real nz=(ccx[2]-midx[2])/dist;

                        Real surfnorm[AMREX_SPACEDIM];

                        surfnorm[0]=(inner_radflag)?-nx:nx;
                        surfnorm[1]=(inner_radflag)?-ny:ny;
                        surfnorm[2]=(inner_radflag)?-nz:nz;

                        Real surf_loc[AMREX_SPACEDIM];

                        //ideally use levelset gradients
                        surf_loc[0]=ccx[0]-(inner_radflag?rad1dist:rad2dist)*surfnorm[0];
                        surf_loc[1]=ccx[1]-(inner_radflag?rad1dist:rad2dist)*surfnorm[1];
                        surf_loc[2]=ccx[2]-(inner_radflag?rad1dist:rad2dist)*surfnorm[2];

                        //amrex::Print()<<"distance, innerrad, outerrad:"<<euclid_dist(surf_loc,midx)<<"\t"
                        //<<inner_radflag<<"\t"<<rad1dist<<"\t"<<rad2dist<<"\n";

                        if(iter==0)
                        {
                            phi_arr(i,j,k)=(inner_radflag)?dircval_rad1:dircval_rad2;
                        }
                        else
                        {
                            //add spacedim ifdefs
                            int nsupport_pts=0;
                            for(int kk=-1;kk<=1;kk++)
                            {
                                for(int jj=-1;jj<=1;jj++)
                                {
                                    for(int ii=-1;ii<=1;ii++)
                                    {
                                        if(!(ii==0 && jj==0 && kk==0))
                                        {
                                            //amrex::Print()<<"ls,ii,jj,kk:"<<ls_arr(i+ii,j+jj,k+kk)<<"\n";
                                            if(ls_arr(i+ii,j+jj,k+kk)==1)
                                            {
                                                stenc_i[nsupport_pts]=ii;
                                                stenc_j[nsupport_pts]=jj;
                                                stenc_k[nsupport_pts]=kk;
                                                nsupport_pts++;
                                            }
                                        }
                                    }
                                }
                            }

                            //amrex::Print()<<"nsupp pts,dx,eps:"<<nsupport_pts<<"\t"<<dx[0]<<"\t"<<eps<<"\n";
                            //amrex::Print()<<"i,j,k:"<<i<<"\t"<<j<<"\t"<<k<<"\n";
                            for(int m=0;m<nsupport_pts;m++)
                            {
                                //amrex::Print()<<"supp pts:"<<stenc_i[m]<<"\t"
                                //<<stenc_j[m]<<"\t"<<stenc_k[m]<<"\n";
                            }

                            for(int m=0;m<nsupport_pts;m++)
                            {
                                Real xm[AMREX_SPACEDIM];

                                xm[0]=prblo[0]+(i+stenc_i[m]+0.5)*dx[0];
                                xm[1]=prblo[1]+(j+stenc_j[m]+0.5)*dx[1];
                                xm[2]=prblo[2]+(k+stenc_k[m]+0.5)*dx[2];

                                bvec[m]=phi_arr(i+stenc_i[m],j+stenc_j[m],k+stenc_k[m]);

                                for(int n=0;n<nsupport_pts;n++)
                                {
                                    Real xn[AMREX_SPACEDIM];

                                    xn[0]=prblo[0]+(i+stenc_i[n]+0.5)*dx[0];
                                    xn[1]=prblo[1]+(j+stenc_j[n]+0.5)*dx[1];
                                    xn[2]=prblo[2]+(k+stenc_k[n]+0.5)*dx[2];

                                    Real r=euclid_dist(xn,xm);

                                    //amrex::Print()<<"m,n,dist:"<<m<<"\t"<<n<<"\t"<<r<<"\n";

                                    mat[m][n]=rbf(r,eps);
                                }
                            }

                            for(int m=0;m<nsupport_pts;m++)
                            {
                                Real xm[AMREX_SPACEDIM];

                                xm[0]=prblo[0]+(i+stenc_i[m]+0.5)*dx[0];
                                xm[1]=prblo[1]+(j+stenc_j[m]+0.5)*dx[1];
                                xm[2]=prblo[2]+(k+stenc_k[m]+0.5)*dx[2];

                                Real r=euclid_dist(xm,surf_loc);
                                //amrex::Print()<<"m,dist:"<<m<<"\t"<<r<<"\n";
                                mat[nsupport_pts][m]=rbf(r,eps);
                                mat[m][nsupport_pts]=mat[nsupport_pts][m];
                            }
                            mat[nsupport_pts][nsupport_pts]=1.0;

                            bvec[nsupport_pts]=(inner_radflag)?dircval_rad1:dircval_rad2;
                            int neq=nsupport_pts+1;

                            //printmat("matrix:", mat, neq);
                            //printvec("bvec:", bvec, neq);

                            bool gmres_succ=performgmres(mat,bvec,soln0,wts,neq,
                                                         0,3,neq,1e-10,0);

                            Real r=euclid_dist(ccx,surf_loc);
                            phi_arr(i,j,k)=wts[nsupport_pts]*rbf(r,eps);
                            for(int m=0;m<nsupport_pts;m++)
                            {
                                Real xm[AMREX_SPACEDIM];

                                xm[0]=prblo[0]+(i+stenc_i[m]+0.5)*dx[0];
                                xm[1]=prblo[1]+(j+stenc_j[m]+0.5)*dx[1];
                                xm[2]=prblo[2]+(k+stenc_k[m]+0.5)*dx[2];
                                Real r=euclid_dist(ccx,xm);
                                phi_arr(i,j,k)+=wts[m]*rbf(r,eps);
                            }
                        }
                        //printvec("wts:", wts, neq);
                        //amrex::Print()<<"phi:"<<"\t"<<phi_arr(i,j,k)<<"\n";
                    }
                });

            }

            for (MFIter mfi(phi); mfi.isValid(); ++mfi)
            {
                const auto dx = geom.CellSizeArray();
                const Box& bx = mfi.tilebox();
                const Box& gbx = amrex::grow(bx, 1);

                Array4<Real> phi_arr = phi.array(mfi);
                Array4<Real> rhs_arr = rhs.array(mfi);
                Array4<Real> ls_arr = levset.array(mfi);
                Array4<Real> bcoef_arr = cell_bcoef.array(mfi);
                Array4<Real> acoef_arr = acoef.array(mfi);
                amrex::Real mindx=std::min(dx[0],std::min(dx[1],dx[2]));
                amrex::Real mindx2=mindx*mindx;

                amrex::ParallelFor(gbx,[=] AMREX_GPU_DEVICE (int i, int j, int k) {

                    if(ls_arr(i,j,k)==1.0)
                    {
                        //rhs_arr(i,j,k) =1.0;
                        rhs_arr(i,j,k) =sourceterm;
                        acoef_arr(i,j,k)=0.0;
                        bcoef_arr(i,j,k)=1.0;
                    }

                    if(ls_arr(i,j,k)>0.0 && ls_arr(i,j,k)<1.0)
                    {
                        bcoef_arr(i,j,k)=1.0;
                        acoef_arr(i,j,k)=penalty/mindx2;
                        rhs_arr(i,j,k)=phi_arr(i,j,k)*penalty/mindx2;
                    }
                    if(ls_arr(i,j,k)==0.0)
                    {
                        rhs_arr(i,j,k)=0.0;
                        bcoef_arr(i,j,k)=1.0;
                        acoef_arr(i,j,k)=penalty/mindx2;
                    }
                });
            }

            mlabec.setACoeffs(0, acoef);
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef), 
                                              cell_bcoef, geom, true);
            mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

            mlmg.solve({&soln}, {&rhs}, tol_rel, tol_abs);
            amrex::MultiFab::Subtract(err,soln, 0, 0, 1 ,0);
            Real errnorm=err.norm2();
            if(iter==0)
            {
                errnorm_1st=errnorm;
            }
            rel_errnorm=errnorm/errnorm_1st;

            amrex::Print()<<"iter, abs errnorm, rel errnorm:"
            <<iter<<"\t"<<errnorm<<"\t"<<rel_errnorm<<"\n";
            amrex::MultiFab::Copy(phi,soln, 0, 0, 1 ,0);
        }

        // store plotfile variables; q and phi
        MultiFab plotfile_mf(grids, dmap, 2, 0);
        MultiFab::Copy(plotfile_mf, levset,0,0,1,0);
        MultiFab::Copy(plotfile_mf, phi,0,1,1,0);

        WriteSingleLevelPlotfile("plt", plotfile_mf, {"levset","potential"}, geom, 0.0, 0);
    }

    amrex::Finalize();
}
