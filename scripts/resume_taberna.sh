#!/usr/bin/env bash
# Reconfigure + build ONLY taberna (build-llvm). Qt and VTK are already installed
# in build-qt/install and build-vtk/install. Clean reconfigure so flag changes
# (PIC, force-includes) take effect.
set -euo pipefail
ROOT=/home/forrest/taberna; cd "$ROOT"
TC="$ROOT/cmake/clang-lld.toolchain.cmake"
LIBCA=/usr/lib/llvm-23/lib/libllvmlibc.a
LLVM_FLAGS=(-DCMAKE_TOOLCHAIN_FILE="$TC" -DTABERNA_LIBCXX=ON
            -DTABERNA_LLVMLIBC=ON -DTABERNA_LLVMLIBC_LIB="$LIBCA"
            -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DCMAKE_BUILD_TYPE=Release)
QT_PREFIX="$ROOT/build-qt/install"

bash "$ROOT/scripts/apply_patches.sh"

echo "=== taberna (clean reconfigure: PIC for shared-lib TLS) ==="
rm -rf build-llvm
cmake -G Ninja "${LLVM_FLAGS[@]}" -DTABERNA_GUI=ON \
  -DCMAKE_PREFIX_PATH="$QT_PREFIX;$ROOT/build-vtk/install" \
  -DVTK_DIR="$ROOT/build-vtk/install/lib/cmake/vtk-9.6" \
  -S . -B build-llvm
ninja -C build-llvm
echo "=== FULL LLVM-STACK BUILD DONE ==="
