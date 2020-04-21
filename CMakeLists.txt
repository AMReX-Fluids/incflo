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
# Import AMReX
#
# If AMREX_HOME is given, AMReX will be imported as a sub-project to be built together with incflo
# If AMREX_HOME is not given, then we will look for a valid AMReX library installation
#
set(AMREX_HOME "" CACHE PATH "Path to AMReX top-level source directory")

if (AMREX_HOME)

  message(STATUS "AMReX top-level source directory path: ${AMREX_HOME}")
  set(AMREX_BINARY_DIR "${incflo_BINARY_DIR}/amrex" CACHE INTERNAL "Path to AMReX top-level binary directory")

  # Add AMReX as a sub-project: incflo inherits all AMReX options
  add_subdirectory(${AMREX_HOME} ${AMREX_BINARY_DIR})

  message(STATUS "AMReX binary directory: ${AMREX_BINARY_DIR}")

  # Add AMReX module path
  list(APPEND CMAKE_MODULE_PATH ${AMREX_HOME}/Tools/CMake)

else()

   #
   # Options -- incflo uses XSDK!!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   #
   option( ENABLE_EB   "Build EB Code" OFF )
   option( ENABLE_MPI  "Enable MPI"  ON)
   option( ENABLE_OMP  "Enable OpenMP" OFF)
   option( ENABLE_CUDA "Enable GPU support via CUDA" OFF )

   #
   # Find amrex
   #
   set(AMREX_MINIMUM_VERSION 20.03 CACHE INTERNAL "Minimum required AMReX version")
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

#
# Add AMReX "build info" capabilities
#
include(AMReXBuildInfo)
generate_buildinfo(incflo ${CMAKE_SOURCE_DIR})


#
# Add cuda support >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#

add_subdirectory(src)
