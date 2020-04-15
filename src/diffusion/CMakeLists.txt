target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   DiffusionScalarOp.cpp
   DiffusionScalarOp.H
   DiffusionTensorOp.cpp
   DiffusionTensorOp.H
   incflo_diffusion.cpp
   )
