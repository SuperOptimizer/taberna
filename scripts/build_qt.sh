#!/usr/bin/env bash
# Minimal Qt 6 from source with the full-LLVM C++ runtime (clang/lld + libc++ +
# libunwind + compiler-rt), same "only the modules we need" approach as VTK. Builds
# ONLY qtbase + qtshadertools + qtdeclarative (Quick/QML/Controls) — not the other
# ~40 Qt modules. Installs to build-qt/install (gitignored). Sources cloned shallow
# into build-qt/src (build artifacts, not tracked submodules).
set -euo pipefail
ROOT=/home/forrest/taberna
cd "$ROOT"
QTVER=v6.10.2                       # match the apt Qt we validated against
PREFIX="$ROOT/build-qt/install"
TC="$ROOT/cmake/clang-lld.toolchain.cmake"
SRC="$ROOT/build-qt/src"
mkdir -p "$SRC"

clone() {  # shallow-clone one Qt module if absent
  [ -e "$SRC/$1/CMakeLists.txt" ] || \
    git clone --depth 1 --branch "$QTVER" "https://github.com/qt/$1.git" "$SRC/$1"
}
clone qtbase
clone qtshadertools
clone qtdeclarative

common=( -G Ninja
  -DCMAKE_TOOLCHAIN_FILE="$TC" -DTABERNA_LIBCXX=ON
  -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_PREFIX_PATH="$PREFIX"
  -DQT_BUILD_TESTS=OFF -DQT_BUILD_EXAMPLES=OFF )

for m in qtbase qtshadertools qtdeclarative; do
  echo "=== configuring $m (libc++) ==="
  cmake "${common[@]}" -S "$SRC/$m" -B "build-qt/$m"
  echo "=== building + installing $m ==="
  ninja -C "build-qt/$m" install
done
echo "=== QT (libc++) MINIMAL BUILD DONE -> $PREFIX ==="
ls "$PREFIX/lib/cmake/" 2>/dev/null | grep -iE "Qt6Quick|Qt6Qml|Qt6Gui" || true
