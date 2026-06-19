# Minimal, dependency-free baseline JPEG reader/writer (own DCT/Huffman tables).
# Consume via add_subdirectory():
#   add_subdirectory(jpg)
#   target_link_libraries(<your_target> jpg)
# which exposes jpg.h on the include path. Standalone (`-DJPG_STANDALONE=ON`)
# also builds the round-trip test.
cmake_minimum_required(VERSION 3.16)
project(jpg C)

option(JPG_STANDALONE "build jpg tests" OFF)
if(CMAKE_PROJECT_NAME STREQUAL "jpg")
  set(JPG_STANDALONE ON)
endif()

add_library(jpg STATIC jpg.c)
target_include_directories(jpg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(jpg PUBLIC m)             # DCT uses libm
set_target_properties(jpg PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(JPG_STANDALONE)
  enable_testing()
  add_executable(jpg_test jpg_test.c)
  target_link_libraries(jpg_test jpg)
  add_test(NAME jpg_roundtrip COMMAND jpg_test)
endif()
