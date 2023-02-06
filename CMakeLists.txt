cmake_minimum_required(VERSION 3.14)

project( incflo
   #VERSION 1.0
   HOMEPAGE_URL "https://amrex-codes.github.io/incflo/"
   DESCRIPTION "An adaptive mesh, variable-density incompressible Navier-Stokes solver with no subcycling in time"
   LANGUAGES C CXX)

#
# Set default build type to Release
#
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()


#
# Options -- incflo uses XSDK!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
set(INCFLO_DIM_VALUES 2 3)
set(INCFLO_DIM 3 CACHE STRING "Dimension of the incflo build. Valid values: <2,3>")
set_property(CACHE INCFLO_DIM PROPERTY STRINGS ${INCFLO_DIM_VALUES})
if(NOT INCFLO_DIM IN_LIST INCFLO_DIM_VALUES)
    message(FATAL_ERROR "\nINCFLO_DIM setto ${INCFLO_DIM}. Must be one of ${INCFLO_DIM_VALUES}. ")
endif()


# We want to offer the user some core options, i.e. MPI, OpenMP, CUDA, HYPRE
# and floating-point exceptions, which require related AMReX components to
# be enabled in the installation.
# In STANDALONE mode, find_package will make sure these options are available
# IN SUPERBUILD mode, These options will supersede AMReX options with the same name
option( INCFLO_OMP    "Enable OpenMP" NO )
option( INCFLO_MPI    "Enable MPI"   YES )
option( INCFLO_CUDA   "Enable CUDA"   NO )
option( INCFLO_EB     "Build Embedded Boundary support" NO )
option( INCFLO_HYPRE  "Enable HYPRE"  NO )
option( INCFLO_FPE    "Enable Floating Point Exceptions checks" NO )

# CSG support
option( INCFLO_CSG "Build with CSG support" OFF)

if (INCFLO_CUDA)
    enable_language(CUDA)
endif ()

#
# Import AMReX
#
# If AMREX_HOME is given, AMReX will be imported as a sub-project to be built together with incflo (SUPERBUILD)
# If AMREX_HOME is not given, then we will look for a valid AMReX library installation (STANDALONE)
#
set(AMREX_HOME "" CACHE PATH "Path to AMReX top-level source directory")

set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "")

if (AMREX_HOME)  # SUPERBUILD MODE

  message(STATUS "AMReX top-level source directory path: ${AMREX_HOME}")
  set(AMREX_BINARY_DIR "${incflo_BINARY_DIR}/amrex" CACHE INTERNAL "Path to AMReX top-level binary directory")

  # Add AMReX module path
  list(APPEND CMAKE_MODULE_PATH ${AMREX_HOME}/Tools/CMake)

  # Set required settings for AMReX
  set(USE_XSDK_DEFAULTS       ON            CACHE INTERNAL "")
  set(XSDK_ENABLE_Fortran     OFF           CACHE INTERNAL "")
  set(XSDK_PRECISION          "DOUBLE"      CACHE INTERNAL "")
  set(AMReX_SPACEDIM          ${INCFLO_DIM} CACHE INTERNAL "")
  set(AMReX_MPI               ${INCFLO_MPI} CACHE INTERNAL "")
  set(AMReX_OMP               ${INCFLO_OMP} CACHE INTERNAL "")
  if (INCFLO_CUDA)
    set(AMReX_GPU_BACKEND     CUDA          CACHE INTERNAL "")
  endif ()
  set(AMReX_EB                ${INCFLO_EB}  CACHE INTERNAL "")
  set(AMREX_FORTRAN           OFF           CACHE INTERNAL "")
  set(AMReX_LINEAR_SOLVERS    ON            CACHE INTERNAL "")
  set(AMReX_BUILD_TUTORIALS   OFF           CACHE INTERNAL "")

  # If CUDA is required, enable the language BEFORE adding the AMReX directory
  # Since AMReX_SetupCUDA has an include guard, it will be included only once here.
  # The reason for enabling CUDA before adding the AMReX subdirectory is that
  # the top-most directory needs to setup the CUDA language before a CUDA-enabled target
  # from a sub-project is included via add_subdirectory.
  # IMPORTANT: if you don't do this, AMReX will perform this step in a sub-scope and therefore
  # it will not setup CUDA here!
  if(INCFLO_CUDA)
     include(AMReX_SetupCUDA)
  endif ()

  # Add AMReX as a sub-project: incflo inherits all AMReX options
  add_subdirectory(${AMREX_HOME} ${AMREX_BINARY_DIR})

  message(STATUS "AMReX binary directory: ${AMREX_BINARY_DIR}")

else ()  # STANDALONE MODE

   # Find amrex
   set(AMREX_MINIMUM_VERSION 20.11 CACHE INTERNAL "Minimum required AMReX version")
   set(AMREX_REQUIRED_COMPONENTS ${INCFLO_DIM}D FORTRAN DOUBLE LSOLVERS)

   if (INCFLO_EB)
      list(APPEND AMREX_REQUIRED_COMPONENTS EB)
   endif ()
   if (INCFLO_MPI)
      list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
   endif ()
   if (INCFLO_OMP)
      list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
   endif ()
   if (INCFLO_CUDA)
      list(APPEND AMREX_REQUIRED_COMPONENTS CUDA)
   endif ()

   find_package(AMReX ${AMREX_MINIMUM_VERSION} CONFIG
      REQUIRED ${AMREX_REQUIRED_COMPONENTS}
      )
   message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")

endif ()

set(AMREX_HYDRO_HOME "" CACHE PATH "Path to AMReX-Hydro top-level source directory")

if (AMREX_HYDRO_HOME)
   message(STATUD "AMReX-Hydro top-level source directory path: ${AMREX_HYDRO_HOME}")
   set(AMREX_HYDRO_BINARY_DIR "${incflo_BINARY_DIR}/AMReX-Hydro" CACHE INTERNAL "Path to AMReX-Hydro top-level binary directory")

   # Set required settings for AMReX-Hydro
   set(HYDRO_SPACEDIM    ${INCFLO_DIM} CACHE STRING "")
   set(HYDRO_MPI         ${INCFLO_MPI} CACHE INTERNAL "")
   set(HYDRO_OMP         ${INCFLO_OMP} CACHE INTERNAL "")
   set(HYDRO_EB          ${INCFLO_EB}  CACHE INTERNAL "")
   if (INCFLO_CUDA)
     set(HYDRO_GPU_BACKEND CUDA CACHE INTERNAL "")
   endif ()

  # Add AMReX-Hydro directory to the build
  add_subdirectory(${AMREX_HYDRO_HOME} ${AMREX_HYDRO_BINARY_DIR})

else ()
   message(FATAL_ERROR "Must provide AMREX_HYDRO_HOME")
endif ()

add_executable(incflo)

#target_link_libraries(incflo PUBLIC AMReX::amrex amrex_hydro)
target_link_libraries(incflo PUBLIC amrex_hydro) #link to amrex through amrex_hydro

#
# If INCFLO_EB=false, we could still link against an external
# EB-enabled AMReX installation since find_package checks for
# presence of a certain component and not for the lack thereof.
# Consequently, the EB-enabled AMReX installation would propagate
# the AMREX_USE_EB definition. Unfortunately, some incflo calls to
# EB routines are guarded by checks on AMREX_USE_EB and are enabled if
# AMReX is EB-enabled. This would generate an error at link time because the
# incflo EB routines get compiled only if INCFLO_EB=yes.
# To workaround this, we use the `-U` preprocessor flag (not sure if
# it's GCC only) to undefine previously defined macros. We do this via
# target_compile_options ecven if it's not designed to handle preprocessor
# definitions.
# The best way forward would be either:
#
#     1) support superbuild only
#     2) guard the EB-related incflo calls to incflo EB-routines via
#        an incflo-specific macro, say INCFLO_USE_EB, dependent on
#        the value of the INCFLO_EB option.
#
if (NOT INCFLO_EB AND AMReX_EB_FOUND)
   target_compile_options(incflo PUBLIC -UAMREX_USE_EB)
endif ()

if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(incflo PUBLIC AMReX::Flags_CXX)
endif ()

if ( INCFLO_FPE )
   target_link_libraries(incflo PUBLIC AMReX::Flags_FPE)
endif ()

#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(incflo ${CMAKE_SOURCE_DIR})


add_subdirectory(src)

set_target_properties(incflo PROPERTIES OUTPUT_NAME "incflo.ex")

#
# Setup CUDA
#
if (INCFLO_CUDA)
   set_cpp_sources_to_cuda_language( incflo )
   set_target_properties( incflo
      PROPERTIES
      CUDA_SEPARABLE_COMPILATION ON  # This add -dc flag
      )
endif ()

if(INCFLO_CSG)
   if(NOT INCFLO_EB)
      message(FATAL_ERROR "CSG option also requires EB.")
   endif()
   if(CSGEB_HOME)
      list(APPEND CMAKE_PREFIX_PATH ${CSGEB_HOME})
   endif()
   find_package(CsgEb)
   target_link_libraries(incflo PUBLIC CsgEb::csg-eb)
endif()
