# clang + lld + full LLVM binutils toolchain for taberna and the vendored VTK.
#
# Prefers an explicitly-versioned upstream LLVM (e.g. the freshly-installed
# llvm-23 binaries) when the suffixed names exist, else falls back to the
# unsuffixed tools on PATH. Pin/override with -DTABERNA_LLVM_SUFFIX=-23 (or "").
#
# Use with:  cmake -DCMAKE_TOOLCHAIN_FILE=cmake/clang-lld.toolchain.cmake ...
# (compiler/linker choice is fixed at first configure — use a fresh build dir
#  when switching a gcc-configured tree to this.)

if(NOT DEFINED TABERNA_LLVM_SUFFIX)
  set(TABERNA_LLVM_SUFFIX "-23")   # try clang-23/llvm-*-23 first, then plain names
endif()

# --- compiler (prefer versioned) ---
find_program(_CLANG   NAMES clang${TABERNA_LLVM_SUFFIX}   clang)
find_program(_CLANGXX NAMES clang++${TABERNA_LLVM_SUFFIX} clang++)
if(_CLANG)
  set(CMAKE_C_COMPILER "${_CLANG}")
endif()
if(_CLANGXX)
  set(CMAKE_CXX_COMPILER "${_CLANGXX}")
endif()

# --- linker: lld via CMake's modern selector (multithreaded by default) ---
set(CMAKE_LINKER_TYPE LLD)
# pin the matching ld.lld if the versioned one exists (else clang's default lld)
find_program(_LD_LLD NAMES ld.lld${TABERNA_LLVM_SUFFIX} ld.lld)
if(_LD_LLD)
  list(APPEND _TABERNA_LD_FLAGS "--ld-path=${_LD_LLD}")
endif()

# --- full LLVM binutils (minimal-gcc: every binary tool is the llvm- one) ---
foreach(_T AR RANLIB NM OBJCOPY OBJDUMP READELF STRIP ADDR2LINE)
  string(TOLOWER "${_T}" _t)
  find_program(_LLVM_${_T} NAMES llvm-${_t}${TABERNA_LLVM_SUFFIX} llvm-${_t})
  if(_LLVM_${_T})
    set(CMAKE_${_T} "${_LLVM_${_T}}")
  endif()
endforeach()

# --- runtime: compiler-rt instead of libgcc. Safe everywhere (builtins are static
# per-object, so this coexists with the gcc-built apt VTK/Qt). We deliberately do
# NOT force --unwindlib=libunwind or -stdlib=libc++ here: the GUI links
# libstdc++/libgcc_s-built Qt+VTK, and a mismatched unwinder/stdlib breaks C++
# exception propagation across those .so boundaries. (Full libc++ requires
# rebuilding Qt+VTK against libc++ — see the toolchain notes.) ---
# Optional FULL-LLVM C++ runtime: libc++ + libc++abi + LLVM libunwind instead of
# libstdc++ + libgcc_s. Only coherent if EVERY C++ lib in the link is libc++ — so
# build Qt and VTK with -DTABERNA_LIBCXX=ON too. Intended for the C++ (GUI/Qt/VTK)
# builds; leave OFF for the pure-C pipeline (no C++ runtime involved there).
if(TABERNA_LIBCXX)
  # libc++'s headers are leaner than libstdc++'s, so code that relied on a name
  # (getenv, std::vector, std::uint32_t, ...) arriving *transitively* fails to
  # compile. Qt and especially VTK have many such latent missing includes. Rather
  # than patch dozens of third-party headers, force-include the common offenders
  # into every C++ TU. Standard headers are include-guarded, so this is a no-op
  # where they're already included; it only fills the gaps.
  set(_TABERNA_FORCE_INC
      "-include cstdlib -include cstring -include cstdint -include cstddef"
      "-include vector -include string -include memory"
      "-include algorithm -include limits -include utility -include functional")
  string(JOIN " " _TABERNA_FORCE_INC ${_TABERNA_FORCE_INC})
  set(CMAKE_CXX_FLAGS_INIT "-stdlib=libc++ ${_TABERNA_FORCE_INC}")
  set(_TABERNA_CXX_RT "-stdlib=libc++" "-unwindlib=libunwind" "-lc++abi")
endif()

# Common link flags for ALL link types: compiler-rt + libc++ runtime + ld path.
string(JOIN " " _TABERNA_COMMON_LINK "--rtlib=compiler-rt" ${_TABERNA_CXX_RT} ${_TABERNA_LD_FLAGS})
set(CMAKE_SHARED_LINKER_FLAGS_INIT "${_TABERNA_COMMON_LINK}")
set(CMAKE_MODULE_LINKER_FLAGS_INIT "${_TABERNA_COMMON_LINK}")

# llvm-libc OVERLAY (apt libllvmlibc-23-dev). EXPERIMENTAL. Its objects use
# local-exec TLS, so the archive CANNOT be linked into shared libraries (Qt/VTK
# .so) — only into final EXECUTABLES. So Qt/VTK build normally (libc++ + glibc) and
# only taberna's own executables get the overlay. (Overlaying Qt/VTK themselves
# would require llvm-libc built as a shared lib, not the static overlay.)
#   -DTABERNA_LLVMLIBC=ON -DTABERNA_LLVMLIBC_LIB=/usr/lib/llvm-23/lib/libllvmlibc.a
set(_TABERNA_EXE_LINK "${_TABERNA_COMMON_LINK}")
if(TABERNA_LLVMLIBC AND TABERNA_LLVMLIBC_LIB)
  string(APPEND _TABERNA_EXE_LINK " ${TABERNA_LLVMLIBC_LIB}")
endif()
set(CMAKE_EXE_LINKER_FLAGS_INIT "${_TABERNA_EXE_LINK}")
