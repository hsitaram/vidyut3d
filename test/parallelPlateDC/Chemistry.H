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

#define N_ID 0
#define I_ID 1
#define NUM_SPECIES (I_ID+1)
#define NUM_ALL_SPECIES (NUM_SPECIES+1)
#define NUM_REACTIONS 1

using namespace amrex;

namespace plasmachem
{
    extern amrex::Vector<std::string> specnames;
    void init();
    void close();
    int find_id(std::string specname);
    
    //FIXME:will this work on AMD gpus?
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real rct_rxnarray[NUM_REACTIONS][NUM_ALL_SPECIES]; 
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real pdt_rxnarray[NUM_REACTIONS][NUM_ALL_SPECIES]; 
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real reaction_elecenergy[NUM_REACTIONS];
    extern AMREX_GPU_DEVICE_MANAGED amrex::Real reaction_gasenergy[NUM_REACTIONS];
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_bg_molwt(amrex::Real specden[NUM_ALL_SPECIES]);
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_charge(int specid);
    
    AMREX_GPU_HOST_DEVICE
    amrex::Real get_molwt(int specid);
    
    AMREX_GPU_HOST_DEVICE 
    void get_reaction_rateconstants(amrex::Real Te,
                                    amrex::Real Tg, amrex::Real Pg,
                                     amrex::Real efield, 
                                     amrex::Real spec[NUM_ALL_SPECIES],
                                    amrex::Real rateconst[NUM_REACTIONS]);
    
    AMREX_GPU_HOST_DEVICE
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg);

    AMREX_GPU_DEVICE
    amrex::Real diffusion_coeff(int specid,
                                amrex::Real Te,
                                amrex::Real efield,
                                amrex::Real Tg, amrex::Real Pg);

    AMREX_GPU_DEVICE
    amrex::Real electron_collision_freq(
                               amrex::Real Te,
                               amrex::Real efield,
                               Real Tg, Real Pg);
}
#endif
