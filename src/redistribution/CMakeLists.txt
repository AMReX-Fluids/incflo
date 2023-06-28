target_include_directories(
   amrex_hydro
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
   )

target_sources(
   amrex_hydro
   PRIVATE
   hydro_redistribution.H
   hydro_redistribution.cpp
   hydro_create_itracker.cpp
   hydro_create_itracker_${HYDRO_SPACEDIM}d_K.H
   hydro_state_redistribute.cpp
   hydro_state_utils.cpp
   )

