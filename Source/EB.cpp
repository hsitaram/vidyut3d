#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB_utils.H>
#include <Vidyut.H>

void Vidyut::init_eb()
{
    // Initialize volume penalization and other complex geometry variables
    std::string geom_kind="all_regular";
    {
        amrex::ParmParse pp("vidyut");
        pp.query("kind_of_geometry",geom_kind);
        
        vp_eps.resize(max_level+1);
        dxmin.resize(max_level+1);
        dxmin2.resize(max_level+1);

        // Find minimum grid spacing to calculate epsilon at each level
        const auto dx = geom[0].CellSizeArray();
        amrex::Real dxmin0 = amrex::min(AMREX_D_DECL(dx[0], dx[1], dx[2]));
        pp.query("vp_epsfact", vp_epsfact);
        pp.query("vp_perm", vp_perm);
        pp.query("vp_iter", vp_iter);

        for(int ilev=0; ilev <= max_level; ilev++){
            dxmin[ilev] = dxmin0/(pow(2.0, ilev));
            dxmin2[ilev] = pow(dxmin[ilev], 2.0);
            vp_eps[ilev] = 1.0 / dxmin[ilev] * vp_epsfact;
        }
    }
    eb_in_domain = true;

    // Generate levelset data at the finest level
    // NOTE: IB implementation assumes refinement factor of 2 for each level
    if(geom_kind == "pins"){

        const int max_pin=2;

        //number of user defined pins
        int num_pin;

        amrex::ParmParse pp("cylinder_pins");
        amrex::Vector<amrex::Array<amrex::Real,AMREX_SPACEDIM>> allcylcent(max_pin);
        amrex::Vector<amrex::Array<amrex::Real,AMREX_SPACEDIM>> allspherecent(max_pin);

        //initalize pins with some dummy values that fall outside of the domain
        const amrex::Real *problo,*probhi;
        amrex::Real maxlen;

        problo=geom[0].ProbLo();
        probhi=geom[0].ProbHi();

        maxlen=amrex::max(AMREX_D_DECL(geom[0].ProbLength(0),geom[0].ProbLength(1), geom[0].ProbLength(2)));

        //setting pins to be way outside the domain initially
        for(int ipin=0;ipin<max_pin;ipin++)
        {
           for(int idim = 0; idim < AMREX_SPACEDIM; idim++) allcylcent[ipin][idim] = problo[idim]-100.0*maxlen;
        }

        //get user defined number of pins (one or two)
        pp.get("num_pin", num_pin);

        // Parsing inputs for pin information, and creating individual cylinder and sphere objects
        amrex::Vector <std::unique_ptr<amrex::EB2::CylinderIF>> impfunc_cylinders(max_pin);
        amrex::Vector <std::unique_ptr<amrex::EB2::SphereIF>> impfunc_spheres(max_pin);
        for(int ipin = 0; ipin < num_pin; ipin++)
        {
            amrex::Array<amrex::Real,AMREX_SPACEDIM> cyl_cent{AMREX_D_DECL(0.0,0.0,0.0)};
            amrex::Array<amrex::Real,AMREX_SPACEDIM> sphere_cent{AMREX_D_DECL(0.0,0.0,0.0)};

            std::string  centstr = "pin_" + std::to_string(ipin) + "_center";
            std::string  rstr = "pin_" + std::to_string(ipin) + "_radius";
            std::string  hstr = "pin_" + std::to_string(ipin) + "_height";
            std::string  dirstr = "pin_" + std::to_string(ipin) + "_dir";
            std::string  signstr = "pin_" + std::to_string(ipin) + "_sign";
            amrex::Vector<amrex::Real> veccent;
            amrex::Real  pinr;
            amrex::Real  pinh;
            int  pindir;
            int  pinsign;
            pp.getarr(centstr.c_str(), veccent,  0, AMREX_SPACEDIM);
            pp.get(rstr.c_str(), pinr);
            pp.get(hstr.c_str(), pinh);
            pp.get(dirstr.c_str(), pindir);
            pp.get(signstr.c_str(), pinsign);
            for(int idir = 0; idir < AMREX_SPACEDIM; idir++)
            {
                cyl_cent[idir] = veccent[idir] ;
                sphere_cent[idir] = veccent[idir];
                if(idir == pindir) sphere_cent[idir] += pinsign * pinh / 2.0;
            }
            allcylcent[ipin] = cyl_cent;
            allspherecent[ipin] = sphere_cent;

            impfunc_cylinders[ipin] = std::unique_ptr<amrex::EB2::CylinderIF>
                                (new amrex::EB2::CylinderIF(pinr, pinh, pindir, allcylcent[ipin], false));
            impfunc_spheres[ipin] = std::unique_ptr<amrex::EB2::SphereIF>
                                (new amrex::EB2::SphereIF(pinr, allspherecent[ipin], false));
        }

        // Make a union of all cylinder and sphere objects to create pin gshop
        auto allpin_IF = amrex::EB2::makeUnion(*impfunc_cylinders[0],*impfunc_cylinders[1],*impfunc_spheres[0],*impfunc_spheres[1]);
        auto gshop = amrex::EB2::makeShop(allpin_IF);

        // Build EB
        amrex::EB2::Build(gshop, geom[max_level], max_level, max_level);
    }
    else
    {
        amrex::Abort("embedded geometry not implemented yet\n");
    }
}
