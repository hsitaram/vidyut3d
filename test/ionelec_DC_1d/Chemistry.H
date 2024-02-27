#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <UnivConstants.H>
#include<VarDefines.H>

/* Elements
0  E
1  He
*/

// Species
#define E_ID 0
#define HE_ID 1
#define HEp_ID 2

#define NUM_ELEMENTS 2
#define NUM_SPECIES 3
#define NUM_IONS 2
#define NUM_REACTIONS 1

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real *awt);
//  MISC
void CKAWT(amrex::Real *awt);
void CKNCF(int *ncf);
void CKSYME_STR(amrex::Vector<std::string> &ename);
void CKSYMS_STR(amrex::Vector<std::string> &kname);
void GET_RMAP(int *_rmap);
void CKINU(const int i, int &nspec, int *ki, int *nu);
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real *x,
            amrex::Real *q_f, amrex::Real *q_r);
void progressRateFR(amrex::Real *q_f, amrex::Real *q_r, amrex::Real *sc,
                    amrex::Real T);
//  SPARSE INFORMATION
void SPARSITY_INFO(int *nJdata, const int *consP, int NCELLS);
void SPARSITY_INFO_SYST(int *nJdata, const int *consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int *nJdata, const int *consP);
void SPARSITY_PREPROC_CSC(int *rowVals, int *colPtrs, const int *consP,
                          int NCELLS);
void SPARSITY_PREPROC_CSR(int *colVals, int *rowPtrs, const int *consP,
                          int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int *colVals, int *rowPtrs, const int *consP,
                               int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int *rowVals, int *colPtrs, int *indx,
                                          const int *consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int *colVals, int *rowPtr,
                                          const int *consP, int base);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int &mm, int &kk, int &ii,
                                                     int &nfit) {
  mm = 2;
  kk = 3;
  ii = 1;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[3] = {
    1822.8884868472639482, // E
    0.2498374807187924,    // HE
    0.2498374807187924,    // HEp
};
#endif
const amrex::Real h_global_imw[3] = {
    1822.8884868472639482, // E
    0.2498374807187924,    // HE
    0.2498374807187924,    // HEp
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[3] = {
    0.000549, // E
    4.002602, // HE
    4.002602, // HEp
};
#endif
const amrex::Real h_global_mw[3] = {
    0.000549, // E
    4.002602, // HE
    4.002602, // HEp
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_imw(amrex::Real *imw_new) {
  imw_new[0] = 1822.8884868472639482; // E
  imw_new[1] = 0.2498374807187924;    // HE
  imw_new[2] = 0.2498374807187924;    // HEp
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real imw(const int n) {
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_mw(amrex::Real *mw_new) {
  mw_new[0] = 0.000549; // E
  mw_new[1] = 4.002602; // HE
  mw_new[2] = 4.002602; // HEp
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real mw(const int n) {
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real &ru, amrex::Real &ruc, amrex::Real &pa) {
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real *species,
                                                   const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];

  // species with no change at a midpoint T
  // species 0: E
  species[0] = +1.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
  // species 1: HE
  species[1] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
  // species 2: HEp
  species[2] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real *species,
                                                   const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];

  // species with no change at a midpoint T
  // species 0: E
  species[0] = +2.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
  // species 1: HE
  species[1] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
  // species 2: HEp
  species[2] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4];
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real *species,
                                                    const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = -7.459783900000000e+02 * invT + 1.423710750000000e+01 -
                 2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.453750000000000e+02 * invT + 1.584651200000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
  } else {
    // species 0: E
    species[0] = -7.459784500000000e+02 * invT + 1.423710750000000e+01 -
                 2.500251500000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.453750000000000e+02 * invT + 1.584651100000000e+00 -
                 2.500000000000000e+00 * tc[0] - 0.000000000000000e+00 * tc[1] -
                 0.000000000000000e+00 * tc[2] - 0.000000000000000e+00 * tc[3] -
                 0.000000000000000e+00 * tc[4];
  }
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real *species,
                                                        const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = -7.45978390e+02 * invT + 1.32371075e+01 -
                 2.50025150e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.45375000e+02 * invT + 5.84651200e-01 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.45375000e+02 * invT + 5.84651200e-01 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
  } else {
    // species 0: E
    species[0] = -7.45978450e+02 * invT + 1.32371075e+01 -
                 2.50025150e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    // species 1: HE
    species[1] = -7.45375000e+02 * invT + 5.84651100e-01 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
    // species 2: HEp
    species[2] = -7.45375000e+02 * invT + 5.84651100e-01 -
                 2.50000000e+00 * tc[0] - 0.00000000e+00 * tc[1] -
                 0.00000000e+00 * tc[2] - 0.00000000e+00 * tc[3] -
                 0.00000000e+00 * tc[4];
  }
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real *species, const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = +1.50025150e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45978390e+02 * invT;
    // species 1: HE
    species[1] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    // species 2: HEp
    species[2] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  } else {
    // species 0: E
    species[0] = +1.50025150e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45978450e+02 * invT;
    // species 1: HE
    species[1] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    // species 2: HEp
    species[2] = +1.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real *species, const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = +2.50025150e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45978390e+02 * invT;
    // species 1: HE
    species[1] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    // species 2: HEp
    species[2] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  } else {
    // species 0: E
    species[0] = +2.50025150e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45978450e+02 * invT;
    // species 1: HE
    species[1] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
    // species 2: HEp
    species[2] = +2.50000000e+00 + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 7.45375000e+02 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real *species, const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 0: E
    species[0] = +2.50025150e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 1.17368560e+01;
    // species 1: HE
    species[1] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 9.15348800e-01;
    // species 2: HEp
    species[2] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 9.15348800e-01;
  } else {
    // species 0: E
    species[0] = +2.50025150e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] - 1.17368560e+01;
    // species 1: HE
    species[1] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 9.15348900e-01;
    // species 2: HEp
    species[2] = +2.50000000e+00 * tc[0] + 0.00000000e+00 * tc[1] +
                 0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
                 0.00000000e+00 * tc[4] + 9.15348900e-01;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real &cpbl) {
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real cpor[3];                           //  temporary storage
  cp_R(cpor, tc);

  // perform dot product
  for (int id = 0; id < 3; ++id) {
    result += x[id] * cpor[id];
  }

  cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real &cpbs) {
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache

  // compute Cp/R at the given temperature

  // species with no change at a midpoint T
  // species 0: E
  result += y[0] *
            (+2.50025150e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            1822.8884868472639482;
  // species 1: HE
  result += y[1] *
            (+2.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            0.2498374807187924;
  // species 2: HEp
  result += y[2] *
            (+2.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            0.2498374807187924;

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real &cvbl) {
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real cvor[3];                           //  temporary storage
  cv_R(cvor, tc);

  // perform dot product
  for (int id = 0; id < 3; ++id) {
    result += x[id] * cvor[id];
  }

  cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real T, const amrex::Real y[], amrex::Real &cvbs) {
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  // compute Cv/R at the given temperature

  // species with no change at a midpoint T
  // species 0: E
  result += y[0] *
            (+1.50025150e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            1822.8884868472639482;
  // species 1: HE
  result += y[1] *
            (+1.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            0.2498374807187924;
  // species 2: HEp
  result += y[2] *
            (+1.50000000e+00 + 0.00000000e+00 * tc[1] + 0.00000000e+00 * tc[2] +
             0.00000000e+00 * tc[3] + 0.00000000e+00 * tc[4]) *
            0.2498374807187924;

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real &hbml) {
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real hml[3];                            //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T
  speciesEnthalpy(hml, tc);

  // perform dot product
  for (int id = 0; id < 3; ++id) {
    result += x[id] * hml[id];
  }

  hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real T, const amrex::Real y[], amrex::Real &hbms) {
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: E
    result += y[0] *
              (+2.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45978390e+02 * invT) *
              1822.8884868472639482;
    // species 1: HE
    result += y[1] *
              (+2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
    // species 2: HEp
    result += y[2] *
              (+2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
  } else {
    // species 0: E
    result += y[0] *
              (+2.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45978450e+02 * invT) *
              1822.8884868472639482;
    // species 1: HE
    result += y[1] *
              (+2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
    // species 2: HEp
    result += y[2] *
              (+2.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
  }

  const amrex::Real RT = 8.31446261815324e+07 * tT; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real &ubml) {
  amrex::Real result = 0;
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real uml[3];                            //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T
  speciesInternalEnergy(uml, tc);

  // perform dot product
  for (int id = 0; id < 3; ++id) {
    result += x[id] * uml[id];
  }

  ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real &ubms) {
  amrex::Real result = 0.0;
  const amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (tT < 1000) {
    // species 0: E
    result += y[0] *
              (+1.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45978390e+02 * invT) *
              1822.8884868472639482;
    // species 1: HE
    result += y[1] *
              (+1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
    // species 2: HEp
    result += y[2] *
              (+1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
  } else {
    // species 0: E
    result += y[0] *
              (+1.50025150e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45978450e+02 * invT) *
              1822.8884868472639482;
    // species 1: HE
    result += y[1] *
              (+1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
    // species 2: HEp
    result += y[2] *
              (+1.50000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3] +
               0.00000000e+00 * tc[4] - 7.45375000e+02 * invT) *
              0.2498374807187924;
  }

  const amrex::Real RT = 8.31446261815324e+07 * tT; // R*T

  ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real x[],
                                                     amrex::Real &sbml) {
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {log(tT), tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[3];                            //  temporary storage
  speciesEntropy(sor, tc);

  // Compute Eq 42
  for (int id = 0; id < 3; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real y[],
                                                     amrex::Real &sbms) {
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {log(tT), tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real sor[3];                            //  temporary storage
  amrex::Real x[3];                              //  need a ytx conversion
  amrex::Real YOW = 0;                           // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (0.000549 * YOW);
  x[1] = y[1] / (4.002602 * YOW);
  x[2] = y[2] / (4.002602 * YOW);
  speciesEntropy(sor, tc);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 3; i++) {
    result += x[i] * (sor[i] - log((x[i] + 1e-100)) - logPratio);
  }
  // Scale by R/W
  sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(const amrex::Real e, const amrex::Real y[], amrex::Real &t,
               int &ierr) {
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real e1, emin, emax, cv, t1, dt;
  CKUBMS(tmin, y, emin);
  CKUBMS(tmax, y, emax);
  if (e < emin) {
    // Linear Extrapolation below tmin
    CKCVBS(tmin, y, cv);
    t = tmin - (emin - e) / cv;
    ierr = 1;
    return;
  }
  if (e > emax) {
    // Linear Extrapolation above tmax
    CKCVBS(tmax, y, cv);
    t = tmax - (emax - e) / cv;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (e - emin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKUBMS(t1, y, e1);
    CKCVBS(t1, y, cv);
    dt = (e - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(const amrex::Real h, const amrex::Real y[], amrex::Real &t,
               int &ierr) {
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real h1, hmin, hmax, cp, t1, dt;
  CKHBMS(tmin, y, hmin);
  CKHBMS(tmax, y, hmax);
  if (h < hmin) {
    // Linear Extrapolation below tmin
    CKCPBS(tmin, y, cp);
    t = tmin - (hmin - h) / cp;
    ierr = 1;
    return;
  }
  if (h > hmax) {
    // Linear Extrapolation above tmax
    CKCPBS(tmax, y, cp);
    t = tmax - (hmax - h) / cp;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (h - hmin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKHBMS(t1, y, h1);
    CKCPBS(t1, y, cp);
    dt = (h - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(const amrex::Real rho,
                                                   const amrex::Real T,
                                                   const amrex::Real x[],
                                                   amrex::Real &P) {
  amrex::Real XW = 0;                      //  To hold mean molecular wt
  XW += x[0] * 0.000549;                   // E
  XW += x[1] * 4.002602;                   // HE
  XW += x[2] * 4.002602;                   // HEp
  P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(const amrex::Real rho,
                                                   const amrex::Real T,
                                                   const amrex::Real y[],
                                                   amrex::Real &P) {
  amrex::Real YOW = 0; //  for computing mean MW

  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }
  P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(const amrex::Real rho,
                                                   const amrex::Real T,
                                                   const amrex::Real c[],
                                                   amrex::Real &P) {
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 0.000549; // E
  W += c[1] * 4.002602; // HE
  W += c[2] * 4.002602; // HEp

  for (int id = 0; id < 3; ++id) {
    sumC += c[id];
  }
  P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real x[],
                                                     amrex::Real &rho) {
  amrex::Real XW = 0;                        //  To hold mean molecular wt
  XW += x[0] * 0.000549;                     // E
  XW += x[1] * 4.002602;                     // HE
  XW += x[2] * 4.002602;                     // HEp
  rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real y[],
                                                     amrex::Real &rho) {
  amrex::Real YOW = 0;

  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }

  rho = P / (8.31446261815324e+07 * T * YOW); // rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real c[],
                                                     amrex::Real &rho) {
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 0.000549; // E
  W += c[1] * 4.002602; // HE
  W += c[2] * 4.002602; // HEp

  for (int id = 0; id < 3; ++id) {
    sumC += c[id];
  }
  rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT(amrex::Real wt[]) {
  get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(const amrex::Real y[],
                                                     amrex::Real &wtm) {
  amrex::Real YOW = 0;

  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(const amrex::Real x[],
                                                     amrex::Real &wtm) {
  amrex::Real XW = 0;    //  see Eq 4 in CK Manual
  XW += x[0] * 0.000549; // E
  XW += x[1] * 4.002602; // HE
  XW += x[2] * 4.002602; // HEp
  wtm = XW;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(const amrex::Real c[],
                                                     amrex::Real &wtm) {
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 0.000549; // E
  W += c[1] * 4.002602; // HE
  W += c[2] * 4.002602; // HEp

  for (int id = 0; id < 3; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against division by zero
  wtm = W / sumC;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(const amrex::Real T,
                                                     amrex::Real cpor[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  cp_R(cpor, tc);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(const amrex::Real T,
                                                     amrex::Real hort[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  speciesEnthalpy(hort, tc);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(const amrex::Real T,
                                                    amrex::Real sor[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {log(tT), tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sor, tc);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(const amrex::Real y[],
                                                    amrex::Real x[]) {
  amrex::Real YOW = 0;

  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 3; i++) {
    x[i] = y[i] * imw(i) * YOWINV;
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real y[],
                                                     amrex::Real c[]) {
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 3; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 3; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 3; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(const amrex::Real rho,
                                                     amrex::Real /*T*/,
                                                     const amrex::Real y[],
                                                     amrex::Real c[]) {

  for (int i = 0; i < 3; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(const amrex::Real x[],
                                                    amrex::Real y[]) {
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 0.000549; // E
  XW += x[1] * 4.002602; // HE
  XW += x[2] * 4.002602; // HEp
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 0.000549 * XWinv;
  y[1] = x[1] * 4.002602 * XWinv;
  y[2] = x[2] * 4.002602 * XWinv;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(const amrex::Real P,
                                                     const amrex::Real T,
                                                     const amrex::Real x[],
                                                     amrex::Real c[]) {
  amrex::Real PORT = P / (8.31446261815324e+07 * T); // P/RT

  // Compute conversion, see Eq 10
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * PORT;
  }
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(const amrex::Real rho,
                                                     const amrex::Real /*T*/,
                                                     const amrex::Real x[],
                                                     amrex::Real c[]) {
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 0.000549; // E
  XW += x[1] * 4.002602; // HE
  XW += x[2] * 4.002602; // HEp
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(const amrex::Real c[],
                                                    amrex::Real x[]) {
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 3; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 3; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(const amrex::Real c[],
                                                    amrex::Real y[]) {
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 0.000549; // E
  CW += c[1] * 4.002602; // HE
  CW += c[2] * 4.002602; // HEp
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 0.000549 * CWinv;
  y[1] = c[1] * 4.002602 * CWinv;
  y[2] = c[2] * 4.002602 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(const amrex::Real T,
                                                     amrex::Real cvml[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  cv_R(cvml, tc);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(const amrex::Real T,
                                                     amrex::Real cpml[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  cp_R(cpml, tc);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(const amrex::Real T,
                                                    amrex::Real uml[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T
  speciesInternalEnergy(uml, tc);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(const amrex::Real T,
                                                    amrex::Real hml[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T
  speciesEnthalpy(hml, tc);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(const amrex::Real T,
                                                    amrex::Real sml[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {log(tT), tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sml, tc);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(const amrex::Real T,
                                                     amrex::Real cvms[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  cv_R(cvms, tc);
  // multiply by R/molecularweight
  cvms[0] *= 1.515633818095350e+11; // E
  cvms[1] *= 2.077264394049980e+07; // HE
  cvms[2] *= 2.077264394049980e+07; // HEp
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(const amrex::Real T,
                                                     amrex::Real cpms[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  cp_R(cpms, tc);
  // multiply by R/molecularweight
  cpms[0] *= 1.515633818095350e+11; // E
  cpms[1] *= 2.077264394049980e+07; // HE
  cpms[2] *= 2.077264394049980e+07; // HEp
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(const amrex::Real T,
                                                    amrex::Real ums[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T

  speciesInternalEnergy(ums, tc);

  for (int i = 0; i < 3; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(const amrex::Real T,
                                                    amrex::Real hms[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {0, tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  amrex::Real RT = 8.31446261815324e+07 * tT;    // R*T

  speciesEnthalpy(hms, tc);

  for (int i = 0; i < 3; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(const amrex::Real T,
                                                    amrex::Real sms[]) {
  amrex::Real tT = T; // temporary temperature
  const amrex::Real tc[5] = {log(tT), tT, tT * tT, tT * tT * tT,
                             tT * tT * tT * tT}; // temperature cache
  speciesEntropy(sms, tc);
  // multiply by R/molecularweight
  sms[0] *= 1.515633818095350e+11; // E
  sms[1] *= 2.077264394049980e+07; // HE
  sms[2] *= 2.077264394049980e+07; // HEp
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specMob(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
    //mobility is scaled by charge
    amrex::Real mob=0.0;
    amrex::Real elecmob=-20.0;
    if(specid==E_ID)
    {
        mob=elecmob;
    }
    if(specid==HEp_ID)
    {
        mob=0.2;
        //mob=0.0;
    }
    if(specid==NUM_SPECIES)
    {
        mob=fivebythree*elecmob;
    }

    return(mob);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
specDiff(const int specid, const amrex::Real Te, const amrex::Real ndens, const amrex::Real emag, const amrex::Real T)
{
    amrex::Real dcoeff=0.0;
    amrex::Real elecdcoeff=100.0;
    if(specid==E_ID)
    {
        dcoeff=elecdcoeff;
    }
    if(specid==HEp_ID)
    {
        dcoeff=0.01;
        //dcoeff=0.0;
    }
    if(specid==NUM_SPECIES)
    {
        dcoeff=fivebythree*elecdcoeff;
    }

    return(dcoeff);
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(amrex::Real *qf, amrex::Real *qr, const amrex::Real *sc,
          const amrex::Real * /*sc_qss*/, const amrex::Real *tc,
          const amrex::Real invT) {

  // reaction 0: HE + E => E + HEp + E
  qf[0] = sc[0] * sc[1];
  qr[0] = 0.0;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 3; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;

  // reaction 0:  HE + E => E + HEp + E
  k_f = 0;
  qf[0] *= k_f;
  qr[0] *= k_f * exp(-(-g_RT[0] + g_RT[1] - g_RT[2])) * (refCinv);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_ener_exch(amrex::Real qf, amrex::Real qr, const amrex::Real * sc, amrex::Real kf, int rxntype, amrex::Real eexci, int
 elidx, amrex::Real * enerExch, amrex::Real Ue, amrex::Real T, const amrex::Real Te)
{
  // (0)='no electron'
  // (1)='excitation/de-excitation'
  // (2)='ionization'
  // (3)='recombination'
  // (4)='attachment'
  // (5)='elastic exchange'
  // (6)='chemi-ioniz'
  // (7)='detachment'

  amrex::Real eV2J = 1.602176620800000e-19;   // eV -> J
  amrex::Real NA = 6.02214085774e23;          // 1/mol
  amrex::Real kB = 1.380649e-23;              // J/K
  amrex::Real ne = sc[E_ID] * NA;             // 1/m3
  amrex::Real mee = Ue / ne;                  // J

  // Molecular weight used to calculate elastic collision energy exchange
  amrex::Real wt[NUM_SPECIES];
  get_mw(wt);

  // Superelastic factor (keeps Te from dropping below gas temperarture)
  amrex::Real sefact = (tanh((Te - T) / 0.5) + 1.0) / 2.0 ;

  // *enerExch = 0.0;
  if(rxntype == 1 || rxntype == 2){
    *enerExch -= (eexci*eV2J)*NA*(qf - qr) * sefact;
  }
  if(rxntype == 3){
    *enerExch -= mee*NA*(qf - qr);
  }
  // if(rxntype[r] == 4){
  //   *enerExch -= mee*NA*(qf[r] - qr[r]);
  // }
  if(rxntype == 5){
    amrex::Real qq = sc[E_ID] * sc[elidx] * kf;
    *enerExch -= 3.0 * kB * (wt[E_ID] / wt[elidx]) * (Te - T) * qq * NA;
  }
  // if(rxntype[r] == 7){
  //   *enerExch += mee*NA*(qf[r] - qr[r]);
  // }

  return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, const amrex::Real* sc, const amrex::Real T, const amrex::Real Te, amrex::Real * enerExch)
{
  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  const amrex::Real invT = 1.0 / tc[1];

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 3; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 3; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // Precalculating values for electron energy exchange evaluation
  amrex::Real ne = sc[E_ID] * 6.02214085774e23;
  amrex::Real Ue = 1.5 * Te * ne * 1.380649e-23;

  {
    // reaction 0:  HE + E => E + HEp + E
    const amrex::Real k_f = 2.5e-12*std::exp(-278508.0/Te)*6.02214076e23;
    const amrex::Real qf = k_f * (sc[0] * sc[1]);
    const amrex::Real qr = 0.0;
    const amrex::Real qdot = qf - qr;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] += qdot;
    int rxntype = 2;
    amrex::Real eexci = 24.0;
    int elidx = 0;
    comp_ener_exch(qf, qr, sc, k_f, rxntype, eexci, elidx, enerExch, Ue, tc[1], Te);
  }

  // Calculate electron energy loss via elastic collisions
  amrex::Real emag = 0.0;
  amrex::Real ndens = 0.0;
  for(int sp=0; sp<NUM_SPECIES; sp++) ndens += sc[sp]*N_A;
  amrex::Real mu = specMob(E_ID,Te,ndens,emag,tc[1]);
  amrex::Real nu = -ECHARGE/ME/mu;
  amrex::Real elec_elastic_coll_term= 3.0/2.0 * K_B * sc[0] * N_A * (Te-tc[1]) * nu * (2.0*ME/(4.0*M_AMU));
  *enerExch -= elec_elastic_coll_term;
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[], const amrex::Real Te, const amrex::Real EN, amrex::Real * ener_exch) {

  // convert to SI
  for (int id = 0; id < 3; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[], 
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[3]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 3; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 3; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(const amrex::Real P,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[], 
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[3]; // temporary storage
  amrex::Real PORT =
      1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real y[],
                                                    amrex::Real wdot[], 
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[3]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 3; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(const amrex::Real rho,
                                                    const amrex::Real T,
                                                    const amrex::Real x[],
                                                    amrex::Real wdot[], 
                                                    const amrex::Real Te,
                                                    amrex::Real * ener_exch) {
  amrex::Real c[3];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 0.000549; // E
  XW += x[1] * 4.002602; // HE
  XW += x[2] * 4.002602; // HEp
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 3; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T, Te, ener_exch);

  // convert to chemkin units
  for (int id = 0; id < 3; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRG(int kcharge[]) {
  kcharge[0] = -1; // E
  kcharge[1] = 0;  // HE
  kcharge[2] = 1;  // HEp
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCHRGMASS(amrex::Real zk[]) {

  int kchrg[3];
  CKCHRG(kchrg);

  for (int id = 0; id < 3; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real *species,
                                                      const amrex::Real *tc) {

  // temperature
  const amrex::Real T = tc[1];

  // species with no change at a midpoint T
  // species 0: E
  species[0] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
  // species 1: HE
  species[1] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
  // species 2: HEp
  species[2] = +0.00000000e+00 + 0.00000000e+00 * tc[1] +
               0.00000000e+00 * tc[2] + 0.00000000e+00 * tc[3];
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(amrex::Real *J, const amrex::Real *sc, const amrex::Real T,
                  const int HP) {

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 16; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[3];
  for (auto &val : wdot) {
    val = 0.0;
  }

  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  amrex::Real invT2 = invT * invT;

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 3; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[3];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[3];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: HE + E => E + HEp + E
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[1];
  k_f = 0;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[0] += q; // E
  wdot[1] -= q; // HE
  wdot[2] += q; // HEp
  // d()/d[E]
  dqdci = +k_f * sc[1];
  J[0] += dqdci; // dwdot[E]/d[E]
  J[1] -= dqdci; // dwdot[HE]/d[E]
  J[2] += dqdci; // dwdot[HEp]/d[E]
  // d()/d[HE]
  dqdci = +k_f * sc[0];
  J[4] += dqdci; // dwdot[E]/d[HE]
  J[5] -= dqdci; // dwdot[HE]/d[HE]
  J[6] += dqdci; // dwdot[HEp]/d[HE]
  // d()/dT
  J[12] += dqdT; // dwdot[E]/dT
  J[13] -= dqdT; // dwdot[HE]/dT
  J[14] += dqdT; // dwdot[HEp]/dT

  amrex::Real c_R[3], dcRdT[3], e_RT[3];
  amrex::Real *eh_RT;
  if (HP == 1) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 3; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[12 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 3; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 3; ++m) {
      dehmixdc += eh_RT[m] * J[k * 4 + m];
    }
    J[k * 4 + 3] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[15] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
#else
  amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(amrex::Real *J, const amrex::Real *sc, const amrex::Real *Tp,
                 const int *HP) {
  amrex::Real c[3];

  for (int k = 0; k < 3; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 3; k++) {
    J[12 + k] *= 1.e-6;
    J[k * 4 + 3] *= 1.e6;
  }
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian(amrex::Real *J,
                                                        const amrex::Real *sc,
                                                        const amrex::Real T,
                                                        const int consP) {

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 16; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[3];
  for (auto &val : wdot) {
    val = 0.0;
  }

  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  amrex::Real invT2 = invT * invT;

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 3; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[3];
  gibbs(g_RT, tc);

  // compute the species enthalpy
  amrex::Real h_RT[3];
  speciesEnthalpy(h_RT, tc);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[3];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: HE + E => E + HEp + E
  // a third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[1];
  k_f = 0;
  dlnkfdT = 0.0;
  // rate of progress
  q = k_f * phi_f;
  dqdT = dlnkfdT * k_f * phi_f;
  // update wdot
  wdot[0] += q; // E
  wdot[1] -= q; // HE
  wdot[2] += q; // HEp
  // d()/d[E]
  dqdci = +k_f * sc[1];
  J[0] += dqdci; // dwdot[E]/d[E]
  J[1] -= dqdci; // dwdot[HE]/d[E]
  J[2] += dqdci; // dwdot[HEp]/d[E]
  // d()/d[HE]
  dqdci = +k_f * sc[0];
  J[4] += dqdci; // dwdot[E]/d[HE]
  J[5] -= dqdci; // dwdot[HE]/d[HE]
  J[6] += dqdci; // dwdot[HEp]/d[HE]
  // d()/dT
  J[12] += dqdT; // dwdot[E]/dT
  J[13] -= dqdT; // dwdot[HE]/dT
  J[14] += dqdT; // dwdot[HEp]/dT

  amrex::Real c_R[3], dcRdT[3], e_RT[3];
  amrex::Real *eh_RT;
  if (consP == 1) {
    cp_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, tc);
    dcvpRdT(dcRdT, tc);
    speciesInternalEnergy(e_RT, tc);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 3; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[12 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 3; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 3; ++m) {
      dehmixdc += eh_RT[m] * J[k * 4 + m];
    }
    J[k * 4 + 3] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[15] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
#else
  amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *J,
                                                    const amrex::Real *sc,
                                                    const amrex::Real *Tp,
                                                    const int *consP) {
  amrex::Real c[3];

  for (int k = 0; k < 3; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 3; k++) {
    J[12 + k] *= 1.e-6;
    J[k * 4 + 3] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetLENIMC(int *LENIMC) {
  *LENIMC = 15;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetLENRMC(int *LENRMC) {
  *LENRMC = 288;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetNO(int *NO) { *NO = 4; }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetKK(int *KK) { *KK = 3; }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetNLITE(int *NLITE) {
  *NLITE = 3;
}

// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetPATM(amrex::Real *PATM) {
  *PATM = 0.1013250000000000E+07;
}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetWT(amrex::Real *WT) {
  WT[0] = 5.48579909E-04;
  WT[1] = 4.00260200E+00;
  WT[2] = 4.00260200E+00;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetEPS(amrex::Real *EPS) {
  EPS[0] = 8.50000000E+02;
  EPS[1] = 1.02000000E+01;
  EPS[2] = 1.02000000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetSIG(amrex::Real *SIG) {
  SIG[0] = 4.25000000E+02;
  SIG[1] = 2.57600000E+00;
  SIG[2] = 2.57600000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetDIP(amrex::Real *DIP) {
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetPOL(amrex::Real *POL) {
  POL[0] = 0.00000000E+00;
  POL[1] = 0.00000000E+00;
  POL[2] = 0.00000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetZROT(amrex::Real *ZROT) {
  ZROT[0] = 1.00000000E+00;
  ZROT[1] = 0.00000000E+00;
  ZROT[2] = 0.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetNLIN(int *NLIN) {
  NLIN[0] = 0;
  NLIN[1] = 0;
  NLIN[2] = 0;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real *COFETA) {
  COFETA[0] = -2.69531463E+01;
  COFETA[1] = -7.88142336E-01;
  COFETA[2] = 3.18038402E-01;
  COFETA[3] = -1.85151697E-02;
  COFETA[4] = -1.39753223E+01;
  COFETA[5] = 1.44641326E+00;
  COFETA[6] = -1.20775609E-01;
  COFETA[7] = 6.06541463E-03;
  COFETA[8] = -1.39753223E+01;
  COFETA[9] = 1.44641326E+00;
  COFETA[10] = -1.20775609E-01;
  COFETA[11] = 6.06541463E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real *COFLAM) {
  COFLAM[0] = 1.12880184E-01;
  COFLAM[1] = -7.88142336E-01;
  COFLAM[2] = 3.18038402E-01;
  COFLAM[3] = -1.85151697E-02;
  COFLAM[4] = 4.19558188E+00;
  COFLAM[5] = 1.44641326E+00;
  COFLAM[6] = -1.20775609E-01;
  COFLAM[7] = 6.06541463E-03;
  COFLAM[8] = 4.19558188E+00;
  COFLAM[9] = 1.44641326E+00;
  COFLAM[10] = -1.20775609E-01;
  COFLAM[11] = 6.06541463E-03;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetCOFD(amrex::Real *COFD) {
  COFD[0] = -1.82216979E+01;
  COFD[1] = 1.68992430E+00;
  COFD[2] = 1.09701574E-01;
  COFD[3] = -8.81888811E-03;
  COFD[4] = -1.70300850E+01;
  COFD[5] = 2.81514717E+00;
  COFD[6] = -1.48815473E-01;
  COFD[7] = 6.35006434E-03;
  COFD[8] = -1.70300850E+01;
  COFD[9] = 2.81514717E+00;
  COFD[10] = -1.48815473E-01;
  COFD[11] = 6.35006434E-03;
  COFD[12] = -1.70300850E+01;
  COFD[13] = 2.81514717E+00;
  COFD[14] = -1.48815473E-01;
  COFD[15] = 6.35006434E-03;
  COFD[16] = -1.19059959E+01;
  COFD[17] = 2.96079381E+00;
  COFD[18] = -1.90702335E-01;
  COFD[19] = 9.23708173E-03;
  COFD[20] = -1.19059959E+01;
  COFD[21] = 2.96079381E+00;
  COFD[22] = -1.90702335E-01;
  COFD[23] = 9.23708173E-03;
  COFD[24] = -1.70300850E+01;
  COFD[25] = 2.81514717E+00;
  COFD[26] = -1.48815473E-01;
  COFD[27] = 6.35006434E-03;
  COFD[28] = -1.19059959E+01;
  COFD[29] = 2.96079381E+00;
  COFD[30] = -1.90702335E-01;
  COFD[31] = 9.23708173E-03;
  COFD[32] = -1.19059959E+01;
  COFD[33] = 2.96079381E+00;
  COFD[34] = -1.90702335E-01;
  COFD[35] = 9.23708173E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void egtransetKTDIF(int *KTDIF) {
  KTDIF[0] = 0;
  KTDIF[1] = 1;
  KTDIF[2] = 2;
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real *COFTD) {
  COFTD[0] = 0.00000000E+00;
  COFTD[1] = 0.00000000E+00;
  COFTD[2] = 0.00000000E+00;
  COFTD[3] = 0.00000000E+00;
  COFTD[4] = 4.35340248E-01;
  COFTD[5] = 1.43259027E-04;
  COFTD[6] = -5.34367083E-08;
  COFTD[7] = 6.20781442E-12;
  COFTD[8] = 4.35340248E-01;
  COFTD[9] = 1.43259027E-04;
  COFTD[10] = -5.34367083E-08;
  COFTD[11] = 6.20781442E-12;
  COFTD[12] = -4.35340248E-01;
  COFTD[13] = -1.43259027E-04;
  COFTD[14] = 5.34367083E-08;
  COFTD[15] = -6.20781442E-12;
  COFTD[16] = 0.00000000E+00;
  COFTD[17] = 0.00000000E+00;
  COFTD[18] = 0.00000000E+00;
  COFTD[19] = 0.00000000E+00;
  COFTD[20] = 0.00000000E+00;
  COFTD[21] = 0.00000000E+00;
  COFTD[22] = 0.00000000E+00;
  COFTD[23] = 0.00000000E+00;
  COFTD[24] = -4.35340248E-01;
  COFTD[25] = -1.43259027E-04;
  COFTD[26] = 5.34367083E-08;
  COFTD[27] = -6.20781442E-12;
  COFTD[28] = 0.00000000E+00;
  COFTD[29] = 0.00000000E+00;
  COFTD[30] = 0.00000000E+00;
  COFTD[31] = 0.00000000E+00;
  COFTD[32] = 0.00000000E+00;
  COFTD[33] = 0.00000000E+00;
  COFTD[34] = 0.00000000E+00;
  COFTD[35] = 0.00000000E+00;
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(amrex::Real *Tci, amrex::Real *ai, amrex::Real *bi,
               amrex::Real *acentric_i) {

  amrex::Real EPS[3];
  amrex::Real SIG[3];
  amrex::Real wt[3];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: E
  Tci[0] = 1.316 * EPS[0];
  ai[0] = (5.55 * avogadro * avogadro * EPS[0] * boltzmann * 1e-24 * SIG[0] *
           SIG[0] * SIG[0]) /
          (wt[0] * wt[0]);
  bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]);
  acentric_i[0] = 0.0;

  // species 1: HE
  Tci[1] = 1.316 * EPS[1];
  ai[1] = (5.55 * avogadro * avogadro * EPS[1] * boltzmann * 1e-24 * SIG[1] *
           SIG[1] * SIG[1]) /
          (wt[1] * wt[1]);
  bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]);
  acentric_i[1] = 0.0;

  // species 2: HEp
  Tci[2] = 1.316 * EPS[2];
  ai[2] = (5.55 * avogadro * avogadro * EPS[2] * boltzmann * 1e-24 * SIG[2] *
           SIG[2] * SIG[2]) /
          (wt[2] * wt[2]);
  bi[2] = 0.855 * avogadro * 1e-24 * SIG[2] * SIG[2] * SIG[2] / (wt[2]);
  acentric_i[2] = 0.0;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(amrex::Real *sqrtOneOverTc, amrex::Real *sqrtAsti,
                   amrex::Real *Bi, amrex::Real *Fomega) {

  // species 0: E
  // Computed from Lennard-Jones
  sqrtOneOverTc[0] = 2.9899408206870e-02;
  sqrtAsti[0] = 7.7623333318521e+12;
  Bi[0] = 7.2051661156636e+10;
  Fomega[0] = 4.8508000000000e-01;

  // species 1: HE
  // Computed from Lennard-Jones
  sqrtOneOverTc[1] = 2.7294300551596e-01;
  sqrtAsti[1] = 5.4994061413731e+04;
  Bi[1] = 2.1989344021800e+00;
  Fomega[1] = 4.8508000000000e-01;

  // species 2: HEp
  // Computed from Lennard-Jones
  sqrtOneOverTc[2] = 2.7294300551596e-01;
  sqrtAsti[2] = 5.4994061413731e+04;
  Bi[2] = 2.1989344021800e+00;
  Fomega[2] = 4.8508000000000e-01;
}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve(amrex::Real * /*A*/, amrex::Real * /*x*/, amrex::Real * /*b*/) {
  amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve_simplified(amrex::Real * /*A*/, amrex::Real * /*x*/,
                    amrex::Real * /*b*/) {
  amrex::Abort(
      "sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
