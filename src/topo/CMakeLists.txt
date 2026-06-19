# taberna topo — native cubical persistent homology (persistence diagrams, later
# Betti matching + topological simplification of the sheetness field).
add_library(taberna_topo STATIC cubical.c)
target_include_directories(taberna_topo PUBLIC ${CMAKE_SOURCE_DIR}/src)
