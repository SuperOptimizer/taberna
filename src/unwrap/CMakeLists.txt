# taberna unwrap — winding field, Archimedean spiral fit, deformation field.
add_library(taberna_unwrap STATIC winding_field.c spiral_fit.c deform.c wmetrics.c)
target_include_directories(taberna_unwrap PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(taberna_unwrap PUBLIC taberna_annotate)
find_library(M_LIB m)
if(M_LIB)
  target_link_libraries(taberna_unwrap PUBLIC ${M_LIB})
endif()

if(OpenMP_C_FOUND)
  target_link_libraries(taberna_unwrap PUBLIC OpenMP::OpenMP_C)
endif()
