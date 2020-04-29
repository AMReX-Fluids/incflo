target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_prob_I.H
   prob_bc.cpp
   prob_bc.H
   prob_init_fluid.cpp
   )
