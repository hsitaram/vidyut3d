#ifndef Vidyut_vardefines_H_
#define Vidyut_vardefines_H_
#include <Chemistry.H>

#define EDN_ID NUM_SPECIES
#define EEN_ID (NUM_SPECIES+1)
#define ETEMP_ID (NUM_SPECIES+2)
#define EDGX_ID (NUM_SPECIES+3)
#define EDGY_ID (NUM_SPECIES+4)
#define EDGZ_ID (NUM_SPECIES+5)
#define POT_ID (NUM_SPECIES+6)
#define EFX_ID (NUM_SPECIES+7)
#define EFY_ID (NUM_SPECIES+8)
#define EFZ_ID (NUM_SPECIES+9)
#define NUM_PLASMAVARS (EFZ_ID-NUM_SPECIES+1)
#define NVAR (NUM_SPECIES+NUM_PLASMAVARS)

#define PERBC 0
#define DIRCBC 1
#define HNEUBC 2
#define IHNEUBC 3
#define ROBINBC 4

#endif
