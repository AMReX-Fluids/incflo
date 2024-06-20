target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   VolumeOfFluid.H
   VolumeOfFluid.cpp
   incflo_vof.cpp
   )
