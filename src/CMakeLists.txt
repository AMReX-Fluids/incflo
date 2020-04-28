target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo.cpp
   incflo_advance.cpp
   incflo_apply_predictor.cpp
   incflo_apply_corrector.cpp
   incflo_compute_dt.cpp
   incflo_compute_forces.cpp
   incflo_tagging.cpp
   incflo_regrid.cpp
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

if (ENABLE_EB)
   add_subdirectory(embedded_boundaries)
endif ()
