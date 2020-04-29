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
   eb_spherecube.cpp
   eb_tuscan.cpp
   eb_twocylinders.cpp
   writeEBsurface.cpp
   eb_if.H
   )
