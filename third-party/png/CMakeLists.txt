# Minimal, dependency-free PNG reader/writer (self-contained DEFLATE/INFLATE).
# Consume via add_subdirectory():
#   add_subdirectory(png)
#   target_link_libraries(<your_target> png)
# which exposes png.h on the include path. Standalone (`-DPNG_STANDALONE=ON`)
# also builds the round-trip test.
#
# NOTE: the target is named `png`. If a sibling dependency (e.g. a bundled
# libpng under VTK) also defines a `png` target, rename this one.
cmake_minimum_required(VERSION 3.16)
project(png C)

option(PNG_STANDALONE "build png tests" OFF)
if(CMAKE_PROJECT_NAME STREQUAL "png")
  set(PNG_STANDALONE ON)
endif()

add_library(png STATIC png.c)
target_include_directories(png PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(png PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(PNG_STANDALONE)
  enable_testing()
  add_executable(png_test png_test.c)
  target_link_libraries(png_test png)
  add_test(NAME png_roundtrip COMMAND png_test)
endif()
