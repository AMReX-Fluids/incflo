target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   incflo_Tracers.cpp
   incflo_PCInit.cpp
   incflo_PCEvolve.cpp
   incflo_PCUtils.cpp
   ParticleData.H
   incflo_PC.H
)

target_compile_definitions(incflo PRIVATE INCFLO_USE_PARTICLES)
