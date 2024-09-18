USE_EB = TRUE
DEBUG = FALSE
USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = FALSE
BL_NO_FORT = TRUE
COMP = llvm

DIM = 3

AMREX_HOME ?= ../../../amrex

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += main.cpp
CEXE_headers += global_defines.H Prob.H

Pdirs := Base Boundary AmrCore
Pdirs += EB
Pdirs += LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

