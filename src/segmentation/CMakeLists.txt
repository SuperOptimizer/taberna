# taberna segmentation — classical (no-ML) volume → graph / sheet primitives.
#   snic.c         — SNIC 3D supervoxels (oversegmentation / graph coarsening)
#   sheet_tensor.c — structure-tensor sheet detection (sheetness + normal field)
#   affinity.c     — signed region-adjacency graph from sheet detection
#   partition.c    — Mutex Watershed signed-graph partition into per-wrap segments
#   ridge.c        — centerline extraction by NMS along the sheet normal
#   stitch.c       — cross-brick label stitching (delayed-merge primitive)
add_library(taberna_segmentation STATIC snic.c sheet_tensor.c ridge.c ced.c trace.c affinity.c partition.c stitch.c)
target_include_directories(taberna_segmentation PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}      # legacy: snic.h included as "snic.h"
  ${CMAKE_SOURCE_DIR}/src)         # new: "segmentation/..", "common/.." style

# libm for the gaussian / eigensolve / sqrt math. No other deps: these operate on
# an in-memory scalar field the caller supplies.
find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(taberna_segmentation PUBLIC ${M_LIB})
endif()

if(OpenMP_C_FOUND)
  target_link_libraries(taberna_segmentation PUBLIC OpenMP::OpenMP_C)
endif()
