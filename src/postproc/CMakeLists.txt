# taberna postproc — classical binary surface-mask cleanup (the Kaggle
# Surface-Detection "win condition": majority/median smoothing, dust removal,
# hole/cavity fill, ball morphology). Depends on taberna_eval for cc_label/topo.
add_library(taberna_postproc STATIC morph.c)
target_include_directories(taberna_postproc PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(taberna_postproc PUBLIC taberna_eval)
