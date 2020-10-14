target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_derive.cpp
   incflo_derive_K.H
   incflo_error.cpp
   )
