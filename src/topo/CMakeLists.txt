# taberna topo — native cubical persistent homology (persistence diagrams, later
# Betti matching + topological simplification of the sheetness field).
add_library(taberna_topo STATIC cubical.c persist0.c)
target_include_directories(taberna_topo PUBLIC ${CMAKE_SOURCE_DIR}/src)

if(OpenMP_C_FOUND)
  target_link_libraries(taberna_topo PUBLIC OpenMP::OpenMP_C)
endif()
