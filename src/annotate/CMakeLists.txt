# taberna annotate — human annotation: umbilicus + winding-annotated point
# collections + link constraints (our own simple text formats).
add_library(taberna_annotate STATIC umbilicus.c point_collection.c)
target_include_directories(taberna_annotate PUBLIC ${CMAKE_SOURCE_DIR}/src)
find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(taberna_annotate PUBLIC ${M_LIB})
endif()
if(OpenMP_C_FOUND)
  target_link_libraries(taberna_annotate PUBLIC OpenMP::OpenMP_C)
endif()
