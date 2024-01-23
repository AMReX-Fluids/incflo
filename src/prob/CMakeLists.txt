target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_prob_I.H
   incflo_prob_usr_I.H
   prob_bc.H
   prob_bc.cpp
   prob_init_fluid.cpp
   prob_init_fluid_usr.cpp
   )
