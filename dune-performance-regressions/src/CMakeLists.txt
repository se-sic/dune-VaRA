add_compile_definitions(GRIDSDIR=\"${CMAKE_CURRENT_SOURCE_DIR}/grids\")

add_executable("dune-performance-regressions" dune-performance-regressions.cc)
add_executable("poisson-test" poisson-test.cc)
add_executable("poisson-alberta" poisson-alberta.cc)
add_executable("poisson-alugrid" poisson-alugrid.cc)
add_executable("poisson-ug-pk-2d" poisson-ug-pk-2d.cc)
add_executable("poisson-yasp-q1-2d" poisson-yasp-q1-2d.cc)
add_executable("poisson-yasp-q2-2d" poisson-yasp-q2-2d.cc)
add_executable("poisson-yasp-q1-3d" poisson-yasp-q1-3d.cc)
add_executable("poisson-yasp-q2-3d" poisson-yasp-q2-3d.cc)

target_link_dune_default_libraries("dune-performance-regressions")
