AMREX_HOME ?= ../../amrex
AMREX_HYDRO_HOME ?= ../../AMReX-Hydro
TOP = ..

EBASE = incflo

DIM = 3

PRECISION  = DOUBLE

TINY_PROFILE = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

USE_HYPRE = FALSE
BL_NO_FORT = TRUE

USE_EB = TRUE

USE_CSG = FALSE

# Set to TRUE to build the tracer-particle module (src/particles, and every
# INCFLO_USE_PARTICLES block); needed for incflo.use_tracer_particles = 1.
USE_PARTICLES = FALSE

include ../src/Make.incflo
