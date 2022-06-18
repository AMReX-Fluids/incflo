target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_projection_bc.cpp
   incflo_apply_nodal_projection.cpp
   )
