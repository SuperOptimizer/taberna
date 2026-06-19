# taberna eval — segmentation / unwrapping metrics (VOI, adapted-Rand, topology).
add_library(taberna_eval STATIC metrics.c topo.c)
target_include_directories(taberna_eval PUBLIC ${CMAKE_SOURCE_DIR}/src)
find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(taberna_eval PUBLIC ${M_LIB})
endif()
