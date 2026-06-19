# taberna viewer — Qt6 + VTK 2x2 volume/MPR frontend.
# Opt-in (the core C pipeline does not depend on Qt/VTK): configure with
#   cmake -S . -B build -DTABERNA_GUI=ON
find_package(Qt6 REQUIRED COMPONENTS Widgets OpenGLWidgets)
find_package(VTK 9 REQUIRED COMPONENTS
  CommonCore CommonDataModel
  RenderingCore RenderingVolume RenderingVolumeOpenGL2 RenderingOpenGL2
  GUISupportQt InteractionStyle)

set(CMAKE_AUTOMOC ON)

add_executable(taberna_viewer
  main.cpp
  main_window.cpp
  slice_view.cpp
  volume_view.cpp
  volume_source.cpp)

target_link_libraries(taberna_viewer PRIVATE
  Qt6::Widgets Qt6::OpenGLWidgets
  ${VTK_LIBRARIES}
  matter_compressor)

# VTK 9 module init (registers the OpenGL2 rendering backend factories)
vtk_module_autoinit(TARGETS taberna_viewer MODULES ${VTK_LIBRARIES})
