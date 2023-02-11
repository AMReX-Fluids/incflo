target_include_directories(incflo PRIVATE ${CMAKE_CURRENT_LIST_DIR})

target_sources(incflo
   PRIVATE
   embedded_boundaries.cpp
   eb_annulus.cpp
   eb_box.cpp
   eb_cylinder.cpp
   eb_cyl_tuscan.cpp
   eb_regular.cpp
   eb_sphere.cpp
   eb_chkptfile.cpp
   eb_if.H)

if (INCFLO_DIM EQUAL 3)
   target_sources(incflo
      PRIVATE
      eb_spherecube.cpp
      eb_tuscan.cpp
      eb_twocylinders.cpp
      writeEBsurface.cpp
)

if (INCFLO_CSG)
   target_sources(incflo PRIVATE eb_csg.cpp)
   target_compile_definitions(incflo PRIVATE CSG_EB)
endif()
