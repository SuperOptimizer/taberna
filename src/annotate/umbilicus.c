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
