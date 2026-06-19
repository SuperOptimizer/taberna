/* tiff_vol.c — see tiff_vol.h. */
#include "io/tiff_vol.h"

#include <stdlib.h>
#include "tiff.h"   // vendored third-party/tiff

#ifdef _OPENMP
#include <omp.h>
#include <unistd.h>
/* OpenMP defaults its team size to the CPU AFFINITY mask, which in many
 * containers/sandboxes is 1 even when many cores are online (nproc=1 but
 * _SC_NPROCESSORS_ONLN=12) — so the parallel loops silently run single-threaded.
 * Every taberna tool links this object (it loads volumes), so this constructor
 * runs at startup and sets the team size to the ONLINE core count instead.
 * Override explicitly with OMP_NUM_THREADS in the environment. */
__attribute__((constructor)) static void taberna_omp_init(void) {
  if (getenv("OMP_NUM_THREADS")) return;          // respect explicit user setting
  long online = sysconf(_SC_NPROCESSORS_ONLN);
  if (online > 1 && online <= 1024) omp_set_num_threads((int)online);
}
#endif

u8 *tiff_load_u8(const char *path, int *nz, int *ny, int *nx) {
  tiff_volume v;
  if (tiff_read_volume(path, &v) != TIFF_OK) return NULL;
  if (v.channels != 1 || v.type != TIFF_U8) { tiff_volume_free(&v); return NULL; }
  if (nz) *nz = (int)v.depth;
  if (ny) *ny = (int)v.height;
  if (nx) *nx = (int)v.width;
  u8 *data = (u8 *)v.data;   // transfer ownership; already z-major u8
  v.data = NULL;
  tiff_volume_free(&v);
  return data;
}

int tiff_save_u8(const char *path, const u8 *data, int nz, int ny, int nx) {
  tiff_volume v = { .width = (unsigned)nx, .height = (unsigned)ny,
                    .depth = (unsigned)nz, .channels = 1, .type = TIFF_U8,
                    .data = (void *)data };
  return tiff_write_volume(path, &v);
}
