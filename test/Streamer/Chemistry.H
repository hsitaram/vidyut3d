#ifndef _CHEMISTRY_H_
#define _CHEMISTRY_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <UnivConstants.H>
#include <ProbParm.H>

#define I_ID 0
#define NUM_SPECIES 1

using namespace amrex;

namespace plasmachem
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
    
    AMREX_GPU_DEVICE AMREX_INLINE
    amrex::Real get_bg_molwt(int i, int j, int k,
                     Array4<Real> const& phi,
                     ProbParm const& prob_parm)
    {
        amrex::Real molwt=40.0*M_AMU;    
        return(molwt);
    }
    
    AMREX_GPU_DEVICE amrex::Real get_charge(int specid);
}
#endif
