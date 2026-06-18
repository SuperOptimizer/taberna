# taberna segmentation — classical (no-ML) volume → graph primitives.
# First member: SNIC 3D supervoxels (snic.c/.h).
add_library(taberna_segmentation STATIC snic.c)
target_include_directories(taberna_segmentation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# libm for the gradient sqr/INFINITY math. No other deps: SNIC operates purely
# on an in-memory scalar field the caller supplies.
find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(taberna_segmentation PUBLIC ${M_LIB})
endif()
