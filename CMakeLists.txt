# taberna — villa-free pipeline to mirror / process / compress Vesuvius scroll
# volumes. Composed of three first-party C submodules under third-party/, wired
# together here so there is ONE shared libs3 across the whole build.
cmake_minimum_required(VERSION 3.20)
project(taberna LANGUAGES C VERSION 0.1.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Frame pointers for profiling: enables the lightweight `perf record --call-graph fp`
# path (frame-pointer unwinding, no costly DWARF unwind). Applied GLOBALLY here, BEFORE
# the third-party subdirs, so the matter-compressor decode hotpath unwinds too.
# Negligible runtime cost on aarch64; off by default (profiling builds set it ON).
option(TABERNA_FRAME_POINTER "Build with frame pointers for perf fp call-graphs" OFF)
if(TABERNA_FRAME_POINTER)
  add_compile_options(-fno-omit-frame-pointer)
endif()

# OpenMP for taberna's compute libs (structure tensor, CED, NSD, morphology, PH):
# the hot loops are per-voxel and embarrassingly parallel. Linked per-target below.
find_package(OpenMP)

# ccache for C and C++, propagated to every subproject. These are directory-scope
# launchers, set (and cached) BEFORE the add_subdirectory()s below, so the libs3 /
# matter-compressor / fysics submodules and the gui inherit them automatically.
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER   "${CCACHE_PROGRAM}" CACHE STRING "ccache (C)"   FORCE)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" CACHE STRING "ccache (C++)" FORCE)
  message(STATUS "ccache enabled for C/C++: ${CCACHE_PROGRAM}")
endif()

# Optional cached ThinLTO for taberna's own targets (NOT VTK, NOT the determinism-
# sensitive submodules — see the apply block below). Off by default: LTO trades
# build time for runtime perf a viewer/research build doesn't need. DTLTO
# (-fthinlto-distributor) is intentionally NOT offered: on a single machine it is
# equivalent to threaded ThinLTO + lld's --thinlto-cache-dir, which we use instead.
option(TABERNA_THINLTO "Cached ThinLTO for taberna's own targets (clang only; release-only)" OFF)
# Fat LTO objects carry BOTH bitcode and machine code, so the .o/.a "gel" with
# non-LTO consumers and plain ar/nm. Costs extra compile time (both codegen paths).
option(TABERNA_FAT_LTO "Emit fat LTO objects when TABERNA_THINLTO is on" ON)

# 1) libs3 — the single shared S3 client (leaf dep). Built first so both
#    matter-compressor and fysics link the SAME `s3` target (no duplicate s3_*
#    symbols). EXCLUDE_FROM_ALL keeps libs3's own tests/examples out of the build.
add_subdirectory(third-party/libs3 EXCLUDE_FROM_ALL)

# 2) matter-compressor — codec + .mca archive. Link the shared s3 instead of its
#    own vendored copy. add_subdirectory => MC_STANDALONE OFF => library only.
set(MC_LIBS3_TARGET s3 CACHE STRING "" FORCE)
add_subdirectory(third-party/matter-compressor EXCLUDE_FROM_ALL)

# 3) fysics — physics kernels + the mca_export pipeline tool. Consumes the
#    `s3` and `matter_compressor` targets defined above.
add_subdirectory(third-party/fysics)

# vendored minimal TIFF reader/writer (now with LZW + multi-page volume read) —
# used by src/io to ingest the Vesuvius surface-detection 320^3 u8 TIFF stacks.
add_subdirectory(third-party/tiff EXCLUDE_FROM_ALL)

# Apply ThinLTO HERE — after the third-party submodules (so it never touches them;
# matter-compressor's byte-determinism forbids LTO-without-PGO) and before our own
# sources, which inherit these directory-scope options. Cache + lld threading give
# the caching/parallelism DTLTO would, on one machine.
if(TABERNA_THINLTO)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    add_compile_options(-flto=thin)
    add_link_options(-flto=thin -Wl,--thinlto-cache-dir=${CMAKE_BINARY_DIR}/thinlto-cache)
    message(STATUS "ThinLTO enabled for taberna targets (cache: ${CMAKE_BINARY_DIR}/thinlto-cache)")
  else()
    message(WARNING "TABERNA_THINLTO needs clang (-flto=thin); ignoring with ${CMAKE_C_COMPILER_ID}")
  endif()
endif()

# Aggressive "ricing" optimization for taberna's OWN targets only — applied HERE,
# after the third-party subdirs, so it NEVER reaches matter-compressor (fast-math
# reorders its DCT and would change decoded bytes; its determinism is sacred). Trades
# IEEE strictness for speed on the FP-heavy winding solve + flatten loops, which the
# profile shows are ~57% of runtime. -march=native isn't otherwise on our targets.
option(TABERNA_FAST "Unsafe-math/native ricing for taberna targets (not third-party)" OFF)
if(TABERNA_FAST)
  add_compile_options(-Ofast -ffast-math -funroll-loops -fno-math-errno
                      -march=native -fno-semantic-interposition)
endif()

# mimalloc — override the system allocator across taberna's OWN targets (not the
# third-party submodules above). Big win on the multithreaded matter-compressor
# decode path (many small concurrent allocs); the winding solver's few huge buffers
# are unaffected. Linked, NOT LD_PRELOAD'd, so it is automatic; --no-as-needed keeps
# the override lib in DT_NEEDED ahead of libc so its malloc/free win. Carried by
# taberna_io (every tool links it) + a constructor in io/mimalloc_tune.c tunes it.
option(TABERNA_MIMALLOC "Override malloc with mimalloc for taberna targets" ON)
add_library(taberna_malloc INTERFACE)
if(TABERNA_MIMALLOC)
  find_library(MIMALLOC_LIB NAMES mimalloc)
  find_path(MIMALLOC_INC NAMES mimalloc.h)
  if(MIMALLOC_LIB AND MIMALLOC_INC)
    target_include_directories(taberna_malloc INTERFACE ${MIMALLOC_INC})
    target_link_libraries(taberna_malloc INTERFACE
      -Wl,--no-as-needed ${MIMALLOC_LIB} -Wl,--as-needed)
    target_compile_definitions(taberna_malloc INTERFACE TABERNA_HAVE_MIMALLOC=1)
    message(STATUS "mimalloc: ${MIMALLOC_LIB}")
  else()
    message(WARNING "TABERNA_MIMALLOC on but mimalloc not found; using system malloc")
  endif()
endif()

# 4) taberna's own sources. segmentation/ turns dense scroll volumes into sheets +
#    the sparse supervoxel graph the un-deformation work is built on; io/ ingests
#    volumes (NRRD now); tools/ are the runnable experiment drivers.
add_subdirectory(src/io)
add_subdirectory(src/annotate)
add_subdirectory(src/segmentation)
add_subdirectory(src/eval)
add_subdirectory(src/topo)
add_subdirectory(src/postproc)
add_subdirectory(src/unwrap)
add_subdirectory(tools)

# 5) optional Qt6+VTK GUI frontend (2x2 MPR + volume render). Opt-in so the core
#    C pipeline builds without Qt/VTK: cmake -DTABERNA_GUI=ON.
option(TABERNA_GUI "Build the Qt6+VTK volume viewer frontend" OFF)
if(TABERNA_GUI)
  enable_language(CXX)            # Qt6/VTK need C++ (core pipeline stays C-only)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_AUTOMOC ON)
  add_subdirectory(gui)
endif()
