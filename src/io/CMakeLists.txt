# taberna io — volume ingest for the experiment phase (NRRD, TIFF stacks; zarr later).
add_library(taberna_io STATIC nrrd.c tiff_vol.c)
target_include_directories(taberna_io PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_include_directories(taberna_io PUBLIC ${CMAKE_SOURCE_DIR}/third-party/tiff)

find_package(ZLIB REQUIRED)
target_link_libraries(taberna_io PUBLIC ZLIB::ZLIB tiff)
