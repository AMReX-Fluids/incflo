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
# We want to offer the user some core options, i.e. MPI, OpenMP, CUDA, HYPRE
# and floating-point exceptions, which require related AMReX components to
# be enabled in the installation.
# In STANDALONE mode, find_package will make sure these options are available
# IN SUPERBUILD mode, These options will supersede AMReX options with the same name
option( ENABLE_OMP    "Enable OpenMP" NO )
option( ENABLE_MPI    "Enable MPI"   YES )
option( ENABLE_HYPRE  "Enable HYPRE"  NO )
option( ENABLE_CUDA   "Enable CUDA"   NO )
option( ENABLE_FPE    "Enable Floating Point Exceptions checks" NO )

if (ENABLE_CUDA)
    enable_language(CUDA)
endif ()

#
# Import AMReX
#
# If AMREX_HOME is given, AMReX will be imported as a sub-project to be built together with incflo (SUPERBUILD)
# If AMREX_HOME is not given, then we will look for a valid AMReX library installation (STANDALONE)
#
set(AMREX_HOME "" CACHE PATH "Path to AMReX top-level source directory")

if (AMREX_HOME)  # SUPERBUILD MODE

  message(STATUS "AMReX top-level source directory path: ${AMREX_HOME}")
  set(AMREX_BINARY_DIR "${incflo_BINARY_DIR}/amrex" CACHE INTERNAL "Path to AMReX top-level binary directory")

  # Add AMReX module path
  list(APPEND CMAKE_MODULE_PATH ${AMREX_HOME}/Tools/CMake)

  # If CUDA is required, enable the language BEFORE adding the AMReX directory
  # Since AMReX_SetupCUDA has an include guard, it will be included only once here.
  # The reason for enabling CUDA before adding the AMReX subdirectory is that
  # the top-most directory needs to setup the CUDA language before a CUDA-enabled target
  # from a sub-project is included via add_subdirectory.
  # IMPORTANT: if you don't do this, AMReX will perform this step in a sub-scope and therefore
  # it will not setup CUDA here!
  if(ENABLE_CUDA)
     include(AMReX_SetupCUDA)
  endif ()

  # Set required settings for AMReX
  set(USE_XSDK_DEFAULTS ON)
  set(DIM 3)
  set(XSDK_PRECISION "DOUBLE" CACHE INTERNAL "Precision:<SINGLE,DOUBLE>" )
  set(ENABLE_LINEAR_SOLVERS ON)
  set(ENABLE_TUTORIALS OFF)

  # Add AMReX as a sub-project: incflo inherits all AMReX options
  add_subdirectory(${AMREX_HOME} ${AMREX_BINARY_DIR})

  message(STATUS "AMReX binary directory: ${AMREX_BINARY_DIR}")

else ()  # STANDALONE MODE

   # Find amrex
   set(AMREX_MINIMUM_VERSION 20.04 CACHE INTERNAL "Minimum required AMReX version")
   set(AMREX_REQUIRED_COMPONENTS 3D DP LSOLVERS)

   if (ENABLE_EB)
      list(APPEND AMREX_REQUIRED_COMPONENTS EB)
   endif ()
   if (ENABLE_MPI)
      list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
   endif ()
   if (ENABLE_OMP)
      list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
   endif ()
   if (ENABLE_CUDA)
      list(APPEND AMREX_REQUIRED_COMPONENTS CUDA)
   endif ()

   find_package(AMReX ${AMREX_MINIMUM_VERSION} CONFIG
      REQUIRED ${AMREX_REQUIRED_COMPONENTS}
      )
   message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")

endif ()

add_executable(incflo)
target_link_libraries(incflo PUBLIC AMReX::amrex)

if ( NOT CMAKE_CXX_FLAGS )
   target_link_libraries(incflo PUBLIC AMReX::Flags_CXX)
endif ()

if ( ENABLE_FPE )
   target_link_libraries(incflo PUBLIC AMReX::Flags_FPE)
endif ()

#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(incflo ${CMAKE_SOURCE_DIR})


add_subdirectory(src)

#
# Setup CUDA
#
if (ENABLE_CUDA)
   set_cpp_sources_to_cuda_language( incflo )
   set_target_properties( incflo
      PROPERTIES
      CUDA_SEPARABLE_COMPILATION ON  # This add -dc flag
      )
endif ()
