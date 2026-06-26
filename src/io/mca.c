/* mca.c — see mca.h. */
#include "io/mca.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "matter_compressor.h"
#include "json.h"

int mca_dims(const char *path, int *nz, int *ny, int *nx, float *quality, int *nlods) {
  int fd = open(path, O_RDONLY);
  if (fd < 0) return -1;
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size <= 0) { close(fd); return -1; }
  size_t len = (size_t)st.st_size;
  void *map = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
  if (map == MAP_FAILED) { close(fd); return -1; }
  mc_reader *r = mc_open((const uint8_t *)map, len);
  if (!r) { munmap(map, len); close(fd); return -1; }
  int x, y, z;
  mc_reader_dims(r, &x, &y, &z);
  if (nx) *nx = x; if (ny) *ny = y; if (nz) *nz = z;
  if (quality) *quality = mc_reader_quality(r);
  if (nlods) *nlods = mc_reader_nlods(r);
  mc_close(r);
  munmap(map, len); close(fd);
  return 0;
}

static mc_archive *open_archive(const char *path) {
  int nx, ny, nz, nlods; float q;
  if (mca_dims(path, &nz, &ny, &nx, &q, &nlods) != 0) return NULL;
  return mc_archive_open_dims(path, nx, ny, nz, q);
}

struct mca_handle {
  mc_archive *a;
  int nz, ny, nx, nlods;
  float q;
};

mca_handle *mca_open(const char *path) {
  int nx, ny, nz, nlods; float q;
  if (mca_dims(path, &nz, &ny, &nx, &q, &nlods) != 0) return NULL;
  mc_archive *a = mc_archive_open_dims(path, nx, ny, nz, q);
  if (!a) return NULL;
  mca_handle *h = (mca_handle *)malloc(sizeof *h);
  if (!h) { mc_archive_close(a); return NULL; }
  h->a = a; h->nz = nz; h->ny = ny; h->nx = nx; h->nlods = nlods; h->q = q;
  return h;
}

void mca_close(mca_handle *h) {
  if (!h) return;
  mc_archive_close(h->a);
  free(h);
}

const char *mca_metadata(const mca_handle *h, size_t *out_len) {
  if (!h) { if (out_len) *out_len = 0; return NULL; }
  size_t len = 0;
  const char *m = mc_archive_metadata(h->a, &len);
  if (out_len) *out_len = m ? len : 0;
  return (m && len) ? m : NULL;
}

int mca_roi_origin(const mca_handle *h, int *oz, int *oy, int *ox) {
  if (oz) *oz = 0; if (oy) *oy = 0; if (ox) *ox = 0;
  size_t len = 0;
  const char *meta = mca_metadata(h, &len);
  if (!meta) return -1;
  json_value *root = json_parse(meta, len);
  if (!root) return -1;
  const json_value *origin = json_obj_get(json_obj_get(root, "roi"), "origin");
  int rc = -1;
  if (origin && origin->type == JSON_ARR && origin->count == 3) {
    if (oz) *oz = json_as_int(json_arr_at(origin, 0), 0);
    if (oy) *oy = json_as_int(json_arr_at(origin, 1), 0);
    if (ox) *ox = json_as_int(json_arr_at(origin, 2), 0);
    rc = 0;
  }
  json_free(root);
  return rc;
}

void mca_handle_dims(const mca_handle *h, int *nz, int *ny, int *nx, float *quality, int *nlods) {
  if (nz) *nz = h->nz; if (ny) *ny = h->ny; if (nx) *nx = h->nx;
  if (quality) *quality = h->q; if (nlods) *nlods = h->nlods;
}

u8 *mca_read(const mca_handle *h, int lod, int z0, int y0, int x0,
             int dz, int dy, int dx) {
  size_t n = (size_t)dz * dy * dx;
  u8 *out = (u8 *)malloc(n);
  if (!out) return NULL;
  size_t sy = (size_t)dx, sz = (size_t)dx * dy;
  mc_archive_read_region(h->a, lod, z0, y0, x0, dz, dy, dx, out, sz, sy, 0);
  return out;
}

u8 *mca_load_region(const char *path, int lod, int z0, int y0, int x0,
                    int dz, int dy, int dx) {
  mc_archive *a = open_archive(path);
  if (!a) return NULL;
  size_t n = (size_t)dz * dy * dx;
  u8 *out = (u8 *)malloc(n);
  if (!out) { mc_archive_close(a); return NULL; }
  size_t sy = (size_t)dx, sz = (size_t)dx * dy;
  mc_archive_read_region(a, lod, z0, y0, x0, dz, dy, dx, out, sz, sy, 0);
  mc_archive_close(a);
  return out;
}

int mca_find_region(const char *path, int lod, int d, float min_frac, uint64_t seed,
                    int *z0, int *y0, int *x0) {
  mc_archive *a = open_archive(path);
  if (!a) return -1;
  mc_box box;
  int got = mc_archive_sample_boxes(a, lod, seed, 1, d, d, d, min_frac, &box);
  mc_archive_close(a);
  if (got < 1) return -1;
  if (z0) *z0 = (int)box.z0; if (y0) *y0 = (int)box.y0; if (x0) *x0 = (int)box.x0;
  return 0;
}
