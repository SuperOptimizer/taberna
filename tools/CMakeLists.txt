# taberna tools — runnable drivers for the experiment phase.
add_executable(sheet_detect sheet_detect.c)
target_link_libraries(sheet_detect PRIVATE taberna_segmentation taberna_io)

add_executable(run_e1 run_e1.c)
target_link_libraries(run_e1 PRIVATE taberna_segmentation taberna_io taberna_eval)

add_executable(extract_cube extract_cube.c)
target_link_libraries(extract_cube PRIVATE taberna_io)

add_executable(eval_seg eval_seg.c)
target_link_libraries(eval_seg PRIVATE taberna_io taberna_eval)

add_executable(run_pipeline run_pipeline.c)
target_link_libraries(run_pipeline PRIVATE
  taberna_io taberna_annotate taberna_segmentation taberna_unwrap)

add_executable(show_inconsistencies show_inconsistencies.c)
target_link_libraries(show_inconsistencies PRIVATE taberna_io taberna_annotate)

add_executable(test_topo test_topo.c)
target_link_libraries(test_topo PRIVATE taberna_eval taberna_postproc)

add_executable(run_surface run_surface.c)
target_link_libraries(run_surface PRIVATE
  taberna_io taberna_segmentation taberna_postproc taberna_eval)

add_executable(surface_sweep surface_sweep.c)
target_link_libraries(surface_sweep PRIVATE
  taberna_io taberna_segmentation taberna_postproc taberna_eval)

add_executable(ridge_clean ridge_clean.c)
target_link_libraries(ridge_clean PRIVATE
  taberna_io taberna_segmentation taberna_postproc taberna_eval)

add_executable(st_sweep st_sweep.c)
target_link_libraries(st_sweep PRIVATE
  taberna_io taberna_segmentation taberna_postproc taberna_eval)

add_executable(ridge_connect ridge_connect.c)
target_link_libraries(ridge_connect PRIVATE
  taberna_io taberna_segmentation taberna_postproc taberna_eval)
