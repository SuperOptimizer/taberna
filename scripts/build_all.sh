#!/usr/bin/env bash
# Single full build of the maximal-LLVM / minimal-GNU stack:
#   clang-23 + lld-23 + llvm binutils + compiler-rt + libc++/libc++abi/libunwind
#   + the apt llvm-libc OVERLAY over system glibc.
# Builds from source against that stack, in order:
#   Qt (qtbase/qtshadertools/qtdeclarative, from official source tarballs)
#   -> VTK (minimal modules) -> taberna (C pipeline + GUI).
# Fail-fast (set -e) with an early smoke test so an unlinkable overlay is caught
# in seconds, not hours. build-*/ is gitignored.
set -euo pipefail
ROOT=/home/forrest/taberna
cd "$ROOT"
TC="$ROOT/cmake/clang-lld.toolchain.cmake"
LIBCA=/usr/lib/llvm-23/lib/libllvmlibc.a       # apt: libllvmlibc-23-dev
[ -f "$LIBCA" ] || { echo "FATAL: $LIBCA missing (apt install libllvmlibc-23-dev)"; exit 1; }

# ---- 0) replay local patches onto vendored submodules (idempotent) ------------
# Submodules check out pristine upstream on a fresh clone; this re-applies our
# libc++/compat fixes so the build is reproducible from a clean tree.
bash "$ROOT/scripts/apply_patches.sh"

LLVM_FLAGS=(-DCMAKE_TOOLCHAIN_FILE="$TC" -DTABERNA_LIBCXX=ON
            -DTABERNA_LLVMLIBC=ON -DTABERNA_LLVMLIBC_LIB="$LIBCA"
            -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DCMAKE_BUILD_TYPE=Release)

# ---- 1) smoke test: does libc++ + the llvm-libc overlay even link + run? -------
echo "=== smoke test: libc++ + llvm-libc overlay ==="
mkdir -p build-smoke
printf '#include <cstdio>\nint main(){ std::puts("ok"); return 0; }\n' > build-smoke/t.cpp
clang++-23 -stdlib=libc++ -unwindlib=libunwind --rtlib=compiler-rt -fuse-ld=lld \
  build-smoke/t.cpp "$LIBCA" -lc++abi -o build-smoke/t \
  && build-smoke/t \
  && echo "SMOKE: full stack links + runs" \
  || { echo "SMOKE FAILED — libc++/overlay link broken; stopping before the big builds"; exit 1; }

# ---- 2) Qt from official source tarballs (minimal: 3 modules) ------------------
QVER=6.10.2
QURL="https://download.qt.io/official_releases/qt/6.10/$QVER/submodules"
SRC="$ROOT/build-qt/src"; mkdir -p "$SRC"
QT_PREFIX="$ROOT/build-qt/install"
fetch_qt() {  # echo the extracted source dir for module $1
  local m=$1 d="$SRC/$m-everywhere-src-$QVER"
  if [ ! -e "$d/CMakeLists.txt" ]; then
    echo "--- fetching $m ---" >&2
    curl -fL --retry 3 -o "$SRC/$m.tar.xz" "$QURL/$m-everywhere-src-$QVER.tar.xz"
    tar -C "$SRC" -xf "$SRC/$m.tar.xz"
  fi
  echo "$d"
}
for m in qtbase qtshadertools qtdeclarative; do
  d="$(fetch_qt "$m")"
  echo "=== Qt: $m ==="
  rm -rf "build-qt/$m"          # fresh configure so CXX_FLAGS_INIT applies (ccache keeps .o fast)
  cmake -G Ninja "${LLVM_FLAGS[@]}" \
    -DCMAKE_INSTALL_PREFIX="$QT_PREFIX" -DCMAKE_PREFIX_PATH="$QT_PREFIX" \
    -DQT_BUILD_TESTS=OFF -DQT_BUILD_EXAMPLES=OFF \
    -S "$d" -B "build-qt/$m"
  ninja -C "build-qt/$m" install
done

# ---- 3) VTK (minimal modules) against our libc++ Qt ---------------------------
echo "=== VTK (libc++ + overlay) ==="
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

# ---- 4) taberna (C pipeline + GUI) against the full stack ----------------------
echo "=== taberna (full stack) ==="
cmake -G Ninja "${LLVM_FLAGS[@]}" -DTABERNA_GUI=ON \
  -DCMAKE_PREFIX_PATH="$QT_PREFIX;$ROOT/build-vtk/install" \
  -DVTK_DIR="$ROOT/build-vtk/install/lib/cmake/vtk-9.6" \
  -S . -B build-llvm
ninja -C build-llvm
echo "=== FULL LLVM-STACK BUILD DONE ==="
