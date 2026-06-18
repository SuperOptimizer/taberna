# taberna — villa-free pipeline to mirror / process / compress Vesuvius scroll
# volumes. Composed of three first-party C submodules under third-party/, wired
# together here so there is ONE shared libs3 across the whole build.
cmake_minimum_required(VERSION 3.20)
project(taberna LANGUAGES C VERSION 0.1.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# 1) libs3 — the single shared S3 client (leaf dep). Built first so both
#    matter-compressor and fysics link the SAME `s3` target (no duplicate s3_*
#    symbols). EXCLUDE_FROM_ALL keeps libs3's own tests/examples out of the build.
add_subdirectory(third-party/libs3 EXCLUDE_FROM_ALL)

# 2) matter-compressor — codec + .mca archive. Link the shared s3 instead of its
#    own vendored copy. add_subdirectory => MC_STANDALONE OFF => library only.
set(MC_LIBS3_TARGET s3 CACHE STRING "" FORCE)
add_subdirectory(third-party/matter-compressor EXCLUDE_FROM_ALL)

# 3) fysics — physics kernels + the mca_export pipeline tool. Consumes the
#    `s3` and `matter_compressor` targets defined above.
add_subdirectory(third-party/fysics)

# 4) taberna's own sources. segmentation/ turns dense scroll volumes into the
#    sparse supervoxel graph the un-deformation work is built on.
add_subdirectory(src/segmentation)
