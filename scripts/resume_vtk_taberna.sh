#!/usr/bin/env bash
# Resume the superbuild from VTK (Qt already installed in build-qt/install).
# Does a CLEAN VTK configure so updated toolchain CXX flags (force-includes) take
# effect, then builds taberna. Qt is left as-is (already built/installed).
set -euo pipefail
ROOT=/home/forrest/taberna; cd "$ROOT"
TC="$ROOT/cmake/clang-lld.toolchain.cmake"
LIBCA=/usr/lib/llvm-23/lib/libllvmlibc.a
LLVM_FLAGS=(-DCMAKE_TOOLCHAIN_FILE="$TC" -DTABERNA_LIBCXX=ON
            -DTABERNA_LLVMLIBC=ON -DTABERNA_LLVMLIBC_LIB="$LIBCA"
            -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
            -DCMAKE_BUILD_TYPE=Release)
QT_PREFIX="$ROOT/build-qt/install"

bash "$ROOT/scripts/apply_patches.sh"

echo "=== VTK (clean reconfigure for new CXX flags) ==="
rm -rf build-vtk
cmake -G Ninja "${LLVM_FLAGS[@]}" \
  -DCMAKE_INSTALL_PREFIX="$ROOT/build-vtk/install" -DCMAKE_PREFIX_PATH="$QT_PREFIX" \
  -DBUILD_SHARED_LIBS=ON -DVTK_BUILD_ALL_MODULES=OFF \
  -DVTK_BUILD_TESTING=OFF -DVTK_BUILD_EXAMPLES=OFF -DVTK_BUILD_DOCUMENTATION=OFF \
  -DVTK_WRAP_PYTHON=OFF -DVTK_WRAP_JAVA=OFF -DVTK_ENABLE_REMOTE_MODULES=OFF \
  -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT -DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
  -DVTK_GROUP_ENABLE_Imaging=DONT_WANT -DVTK_GROUP_ENABLE_Views=DONT_WANT \
  -DVTK_GROUP_ENABLE_Qt=DONT_WANT -DVTK_GROUP_ENABLE_MPI=NO -DVTK_GROUP_ENABLE_Web=NO \
  -DVTK_QT_VERSION=6 \
  -DVTK_MODULE_ENABLE_VTK_GUISupportQtQuick=YES -DVTK_MODULE_ENABLE_VTK_GUISupportQt=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingVolumeOpenGL2=YES -DVTK_MODULE_ENABLE_VTK_RenderingOpenGL2=YES \
  -DVTK_MODULE_ENABLE_VTK_RenderingUI=YES -DVTK_MODULE_ENABLE_VTK_InteractionStyle=YES \
  -S third-party/vtk -B build-vtk
ninja -C build-vtk install

echo "=== taberna (full stack) ==="
cmake -G Ninja "${LLVM_FLAGS[@]}" -DTABERNA_GUI=ON \
  -DCMAKE_PREFIX_PATH="$QT_PREFIX;$ROOT/build-vtk/install" \
  -DVTK_DIR="$ROOT/build-vtk/install/lib/cmake/vtk-9.6" \
  -S . -B build-llvm
ninja -C build-llvm
echo "=== FULL LLVM-STACK BUILD DONE ==="
