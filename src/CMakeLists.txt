target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo.cpp
   incflo_advance.cpp
   incflo_apply_predictor.cpp
   incflo_apply_corrector.cpp
   incflo_compute_dt.cpp
   incflo_compute_forces.cpp
   incflo_correct_small_cells.cpp
   incflo_redistribute.cpp
   incflo_explicit_update.cpp
   incflo_regrid.cpp
   incflo_tagging.cpp
   incflo_update_density.cpp
   incflo_update_tracer.cpp
   incflo_update_velocity.cpp
   incflo_utils.cpp
   main.cpp
   )

add_subdirectory(boundary_conditions)
add_subdirectory(convection)
add_subdirectory(derive)
add_subdirectory(diffusion)
add_subdirectory(prob)
add_subdirectory(projection)
add_subdirectory(rheology)
add_subdirectory(setup)
add_subdirectory(utilities)

if (INCFLO_EB)
   add_subdirectory(embedded_boundaries)
endif ()
