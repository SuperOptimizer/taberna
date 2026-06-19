/* umbilicus.c — see umbilicus.h. */
#include "annotate/umbilicus.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int umbilicus_load(const char *path, umbilicus *u) {
  memset(u, 0, sizeof *u);
  FILE *f = fopen(path, "r");
  if (!f) return -1;
  int cap = 16;
  u->z = malloc(cap * sizeof(f32));
  u->y = malloc(cap * sizeof(f32));
  u->x = malloc(cap * sizeof(f32));
  char line[256];
  while (fgets(line, sizeof line, f)) {
    char *p = line;
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '#' || *p == '\n' || *p == '\0') continue;
    float z, y, x;
    if (sscanf(p, "%f %f %f", &z, &y, &x) != 3) continue;
    if (u->n == cap) {
      cap *= 2;
      u->z = realloc(u->z, cap * sizeof(f32));
      u->y = realloc(u->y, cap * sizeof(f32));
      u->x = realloc(u->x, cap * sizeof(f32));
    }
    u->z[u->n] = z; u->y[u->n] = y; u->x[u->n] = x; u->n++;
  }
  fclose(f);

  // insertion sort by z (control-point counts are tiny)
  for (int i = 1; i < u->n; i++) {
    f32 z = u->z[i], y = u->y[i], x = u->x[i];
    int j = i - 1;
    while (j >= 0 && u->z[j] > z) {
      u->z[j + 1] = u->z[j]; u->y[j + 1] = u->y[j]; u->x[j + 1] = u->x[j]; j--;
    }
    u->z[j + 1] = z; u->y[j + 1] = y; u->x[j + 1] = x;
  }
  return u->n > 0 ? 0 : -2;
}

void umbilicus_free(umbilicus *u) {
  free(u->z); free(u->y); free(u->x);
  memset(u, 0, sizeof *u);
}

int umbilicus_save(const char *path, const umbilicus *u) {
  FILE *f = fopen(path, "w");
  if (!f) return -1;
  fprintf(f, "# z y x  (taberna umbilicus, auto-estimated)\n");
  for (int i = 0; i < u->n; i++)
    fprintf(f, "%.2f %.2f %.2f\n", u->z[i], u->y[i], u->x[i]);
  fclose(f);
  return 0;
}

// Angular variance of the material's outer radial extent seen from (cy,cx) in
// slice z. The scroll's outer boundary is ~circular about the true axis, so the
// center is the point minimizing this variance.
static double umb_sym_score(const u8 *mask, int ny, int nx, int z,
                            double cy, double cx) {
  enum { NA = 72 };
  double rmax[NA];
  int maxr = (int)(0.5 * (ny < nx ? ny : nx));
  size_t nynx = (size_t)ny * nx;
  for (int a = 0; a < NA; a++) {
    double an = a * 2.0 * M_PI / NA, ca = cos(an), sa = sin(an);
    double far = 0;
    for (int r = 1; r < maxr; r++) {
      int yy = (int)(cy + r * sa), xx = (int)(cx + r * ca);
      if (yy < 0 || yy >= ny || xx < 0 || xx >= nx) break;
      if (mask[(size_t)z * nynx + (size_t)yy * nx + xx]) far = r;
    }
    rmax[a] = far;
  }
  double m = 0;
  for (int a = 0; a < NA; a++) m += rmax[a];
  m /= NA;
  double v = 0;
  for (int a = 0; a < NA; a++) { double d = rmax[a] - m; v += d * d; }
  return v / NA;
}

int umbilicus_estimate(const u8 *mask, int nz, int ny, int nx, int nctrl,
                       umbilicus *out) {
  if (nctrl < 2) nctrl = 2;
  memset(out, 0, sizeof *out);
  out->n = nctrl;
  out->z = malloc(nctrl * sizeof(f32));
  out->y = malloc(nctrl * sizeof(f32));
  out->x = malloc(nctrl * sizeof(f32));
  if (!out->z || !out->y || !out->x) { umbilicus_free(out); return -1; }
  size_t nynx = (size_t)ny * nx;

  for (int b = 0; b < nctrl; b++) {
    int z = (int)((b + 0.5) * nz / nctrl);
    if (z >= nz) z = nz - 1;
    // centroid initial guess
    double sy = 0, sx = 0;
    size_t cnt = 0;
    #pragma omp parallel for schedule(static) reduction(+:sy,sx,cnt)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++)
        if (mask[(size_t)z * nynx + (size_t)y * nx + x]) { sy += y; sx += x; cnt++; }
    double cy = cnt ? sy / cnt : ny * 0.5, cx = cnt ? sx / cnt : nx * 0.5;
    // refine: coarse-to-fine local search minimizing angular variance
    double step = 0.1 * (ny < nx ? ny : nx);
    for (int iter = 0; iter < 7; iter++) {
      double best = umb_sym_score(mask, ny, nx, z, cy, cx), bcy = cy, bcx = cx;
      for (int dyi = -1; dyi <= 1; dyi++)
        for (int dxi = -1; dxi <= 1; dxi++) {
          if (!dyi && !dxi) continue;
          double ty = cy + dyi * step, tx = cx + dxi * step;
          double s = umb_sym_score(mask, ny, nx, z, ty, tx);
          if (s < best) { best = s; bcy = ty; bcx = tx; }
        }
      cy = bcy; cx = bcx; step *= 0.5;
    }
    out->z[b] = (f32)z; out->y[b] = (f32)cy; out->x[b] = (f32)cx;
  }

  // the axis is smooth in z: 3-point moving average on (y,x)
  if (nctrl >= 3) {
    f32 *ty = malloc(nctrl * sizeof(f32)), *tx = malloc(nctrl * sizeof(f32));
    for (int i = 0; i < nctrl; i++) {
      int a = i > 0 ? i - 1 : i, c = i < nctrl - 1 ? i + 1 : i;
      ty[i] = (out->y[a] + out->y[i] + out->y[c]) / 3.0f;
      tx[i] = (out->x[a] + out->x[i] + out->x[c]) / 3.0f;
    }
    memcpy(out->y, ty, nctrl * sizeof(f32));
    memcpy(out->x, tx, nctrl * sizeof(f32));
    free(ty); free(tx);
  }
  return 0;
}

void umbilicus_center(const umbilicus *u, f32 z, f32 *cy, f32 *cx) {
  if (u->n == 0) { *cy = 0; *cx = 0; return; }
  if (u->n == 1 || z <= u->z[0]) { *cy = u->y[0]; *cx = u->x[0]; return; }
  if (z >= u->z[u->n - 1]) { *cy = u->y[u->n - 1]; *cx = u->x[u->n - 1]; return; }
  int i = 0;
  while (i < u->n - 1 && u->z[i + 1] < z) i++;
  f32 z0 = u->z[i], z1 = u->z[i + 1];
  f32 t = (z1 > z0) ? (z - z0) / (z1 - z0) : 0.0f;
  *cy = u->y[i] + t * (u->y[i + 1] - u->y[i]);
  *cx = u->x[i] + t * (u->x[i + 1] - u->x[i]);
}

void umbilicus_polar(const umbilicus *u, f32 z, f32 y, f32 x, f32 *theta, f32 *radius) {
  f32 cy, cx;
  umbilicus_center(u, z, &cy, &cx);
  f32 dy = y - cy, dx = x - cx;
  *theta = atan2f(dy, dx);
  *radius = sqrtf(dy * dy + dx * dx);
}
