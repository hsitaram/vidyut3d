#include "Chemistry.H"
const int rmap[9] = {2, 3, 4, 7, 8, 0, 1, 5, 6};

// Returns 0-based map of reaction order
void GET_RMAP(int *_rmap) {
  for (int j = 0; j < 9; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int &nspec, int ki[], int nu[]) {
  const int ns[9] = {4, 4, 3, 3, 3, 4, 4, 3, 3};
  const int kiv[36] = {0, 1, 0, 3, 0, 3, 0, 1, 1, 0, 2, 0, 3, 0, 2, 0, 5, 0,
                       4, 0, 0, 4, 1, 3, 3, 0, 1, 2, 1, 3, 5, 0, 1, 2, 4, 0};
  const int nuv[36] = {-1, -1, 1, 1, -1, -1, 1, 1, -1, 1,  1, 0,
                       -1, 1,  1, 0, -1, 1,  1, 0, -1, -1, 1, 1,
                       -2, 1,  1, 1, -1, -1, 1, 0, -1, -1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 4;
  } else {
    if (i > 9) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 4 + j] + 1;
        nu[j] = nuv[(i - 1) * 4 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(const amrex::Real P, const amrex::Real T, const amrex::Real x[],
            amrex::Real q_f[], amrex::Real q_r[]) {
  amrex::Real c[6]; // temporary storage
  amrex::Real PORT =
      1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 6; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 9; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *q_f, amrex::Real *q_r, amrex::Real *sc,
                    amrex::Real T) {
  const amrex::Real tc[5] = {log(T), T, T * T, T * T * T,
                             T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  // compute the Gibbs free energy
  amrex::Real g_RT[6];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void atomicWeight(amrex::Real *awt) {
  awt[0] = 0.000549; // E
  awt[1] = 4.002602; // He
}

// get atomic weight for all elements
void CKAWT(amrex::Real *awt) { atomicWeight(awt); }

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void CKNCF(int *ncf) {
  int kd = 2;
  // Zero ncf
  for (int id = 0; id < kd * 6; ++id) {
    ncf[id] = 0;
  }

  // E
  ncf[0 * kd + 0] = 1; // E

  // HE
  ncf[1 * kd + 1] = 1; // He

  // HEp
  ncf[2 * kd + 0] = -1; // E
  ncf[2 * kd + 1] = 1;  // He

  // HEm
  ncf[3 * kd + 1] = 1; // He

  // HE2p
  ncf[4 * kd + 0] = -1; // E
  ncf[4 * kd + 1] = 2;  // He

  // HE2m
  ncf[5 * kd + 1] = 2; // He
}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string> &ename) {
  ename.resize(2);
  ename[0] = "E";
  ename[1] = "He";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string> &kname) {
  kname.resize(6);
  kname[0] = "E";
  kname[1] = "HE";
  kname[2] = "HEp";
  kname[3] = "HEm";
  kname[4] = "HE2p";
  kname[5] = "HE2m";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO(int *nJdata, const int *consP, int NCELLS) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 7; k++) {
    for (int l = 0; l < 7; l++) {
      if (Jac[7 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST(int *nJdata, const int *consP, int NCELLS) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 7; k++) {
    for (int l = 0; l < 7; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[7 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED(int *nJdata, const int *consP) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 7; k++) {
    for (int l = 0; l < 7; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[7 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void SPARSITY_PREPROC_CSC(int *rowVals, int *colPtrs, const int *consP,
                          int NCELLS) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 7;
    int offset_col = nc * 7;
    for (int k = 0; k < 7; k++) {
      for (int l = 0; l < 7; l++) {
        if (Jac[7 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void SPARSITY_PREPROC_CSR(int *colVals, int *rowPtrs, const int *consP,
                          int NCELLS, int base) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 7;
      for (int l = 0; l < 7; l++) {
        for (int k = 0; k < 7; k++) {
          if (Jac[7 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 7;
      for (int l = 0; l < 7; l++) {
        for (int k = 0; k < 7; k++) {
          if (Jac[7 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int *colVals, int *rowPtr, const int *consP,
                               int NCELLS, int base) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 7;
      for (int l = 0; l < 7; l++) {
        for (int k = 0; k < 7; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[7 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 7;
      for (int l = 0; l < 7; l++) {
        for (int k = 0; k < 7; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[7 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int *rowVals, int *colPtrs, int *indx,
                                          const int *consP) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 7; k++) {
    for (int l = 0; l < 7; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 7 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[7 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 7 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int *colVals, int *rowPtr,
                                          const int *consP, int base) {
  amrex::GpuArray<amrex::Real, 49> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 6> conc = {0.0};
  for (int n = 0; n < 6; n++) {
    conc[n] = 1.0 / 6.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 7; l++) {
      for (int k = 0; k < 7; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[7 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 7; l++) {
      for (int k = 0; k < 7; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[7 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}