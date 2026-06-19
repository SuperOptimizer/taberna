#!/usr/bin/env bash
# Build a MINIMAL VTK (only the modules taberna's QML viewer needs) from the
# vendored submodule, with clang/lld + ccache. Installs to build-vtk/install.
# Re-run is cheap (ccache + ninja incremental). build-vtk/ is gitignored.
set -euo pipefail
ROOT=/home/forrest/taberna
cd "$ROOT"

# 1) VTK as a shallow submodule on master (latest), if not already present.
if [ ! -e third-party/vtk/CMakeLists.txt ]; then
  echo "=== adding VTK submodule (shallow, master) ==="
  git submodule add --depth 1 https://github.com/Kitware/VTK.git third-party/vtk
  git config -f .gitmodules submodule.third-party/vtk.shallow true
fi

# 2) Configure: disable EVERYTHING, then enable only the few modules we need.
#    GUISupportQtQuick -> QQuickVTKItem (VTK-in-QML); RenderingVolumeOpenGL2 -> GPU
#    volume ray-cast. Dependency closure pulls in CommonCore/RenderingCore/etc.
echo "=== configuring minimal VTK ==="
cmake -G Ninja -S third-party/vtk -B build-vtk \
  -DCMAKE_TOOLCHAIN_FILE="$ROOT/cmake/clang-lld.toolchain.cmake" \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$ROOT/build-vtk/install" \
  -DBUILD_SHARED_LIBS=ON \
  -DVTK_BUILD_ALL_MODULES=OFF \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_BUILD_DOCUMENTATION=OFF \
  -DVTK_WRAP_PYTHON=OFF \
  -DVTK_WRAP_JAVA=OFF \
  -DVTK_ENABLE_REMOTE_MODULES=OFF \
  -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT \
  -DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
  -DVTK_GROUP_ENABLE_Imaging=DONT_WANT \
  -DVTK_GROUP_ENABLE_Views=DONT_WANT \
  -DVTK_GROUP_ENABLE_Qt=DONT_WANT \
  -DVTK_GROUP_ENABLE_MPI=NO \
  -DVTK_GROUP_ENABLE_Web=NO \
  -DVTK_QT_VERSION=6 \
  -DVTK_MODULE_ENABLE_VTK_GUISupportQtQuick=YES \
  -DVTK_MODULE_ENABLE_VTK_GUISupportQt=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingVolumeOpenGL2=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingOpenGL2=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingUI=YES \
  -DVTK_MODULE_ENABLE_VTK_InteractionStyle=YES

# 3) Build + install.
echo "=== building VTK ==="
ninja -C build-vtk install
echo "=== VTK BUILD+INSTALL DONE -> $ROOT/build-vtk/install ==="
ls "$ROOT/build-vtk/install/lib/cmake/" 2>/dev/null
