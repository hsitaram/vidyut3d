#ifndef echemAMR_constants_H_
#define echemAMR_constants_H_
#include <Chemistry.H>

#define N_A 6.02214076e23
#define K_B 1.380649e-23
#define PI 3.14159265358979323846264338327950288
#define P_NTP 101325.0 //pressure at normal pressure and temperature
#define T_NTP 293.15  //temperature at normal pressure and temperature
#define FARADCONST 96485.332
#define GASCONST   8.31446261815324
#define ECHARGE    1.60217662e-19
#define EPS0       8.854187817e-12
#define eV ECHARGE/K_B
#define FBYR FARADCONST/GASCONST

#define EDN_ID NUM_SPECIES
#define EEN_ID (NUM_SPECIES+1)
#define POT_ID (NUM_SPECIES+2)
#define EFX_ID (NUM_SPECIES+3)
#define EFY_ID (NUM_SPECIES+4)
#define EFZ_ID (NUM_SPECIES+5)
#define NUM_PLASMAVARS (EFZ_ID-NUM_SPECIES+1)
#define NVAR (NUM_SPECIES+NUM_PLASMAVARS)

#endif
