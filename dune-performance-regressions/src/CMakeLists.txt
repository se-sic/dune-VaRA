add_executable("dune-performance-regressions" dune-performance-regressions.cc)
add_executable("poisson-test" poisson-test.cc)
target_link_dune_default_libraries("dune-performance-regressions")
