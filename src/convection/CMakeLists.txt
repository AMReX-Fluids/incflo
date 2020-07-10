target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_compute_advection_term.cpp
   incflo_correct_small_cells.cpp
   incflo_MAC_projection.cpp
   incflo_mol_predict_eb.cpp
   incflo_mol_predict.cpp
   incflo_mol_fluxes_eb.cpp
   incflo_mol_fluxes.cpp
   incflo_godunov_predict_${DIM}D.cpp
   incflo_godunov_advection_${DIM}D.cpp
   incflo_godunov_plm.cpp
   incflo_godunov_ppm.cpp
   incflo_godunov_plm.H
   incflo_godunov_ppm.H
   )
