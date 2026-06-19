/* point_collection.c — see point_collection.h. */
#include "annotate/point_collection.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static point_collection *pc_new_collection(pc_set *s, const char *id, int absolute) {
  if (s->ncols == s->cols_cap) {
    s->cols_cap = s->cols_cap ? s->cols_cap * 2 : 8;
    s->cols = realloc(s->cols, (size_t)s->cols_cap * sizeof(point_collection));
  }
  point_collection *c = &s->cols[s->ncols++];
  memset(c, 0, sizeof *c);
  strncpy(c->id, id, sizeof c->id - 1);
  c->absolute = absolute;
  return c;
}

static void pc_add_point(point_collection *c, f32 z, f32 y, f32 x, f32 w) {
  if (c->npts == c->cap) {
    c->cap = c->cap ? c->cap * 2 : 8;
    c->pts = realloc(c->pts, (size_t)c->cap * sizeof(anno_point));
  }
  c->pts[c->npts++] = (anno_point){z, y, x, w};
}

static void pc_add_link(pc_set *s, link_constraint lk) {
  if (s->nlinks == s->links_cap) {
    s->links_cap = s->links_cap ? s->links_cap * 2 : 8;
    s->links = realloc(s->links, (size_t)s->links_cap * sizeof(link_constraint));
  }
  s->links[s->nlinks++] = lk;
}

int pc_load(const char *path, pc_set *s) {
  memset(s, 0, sizeof *s);
  FILE *f = fopen(path, "r");
  if (!f) return -1;
  char line[512];
  point_collection *cur = NULL;
  while (fgets(line, sizeof line, f)) {
    char *p = line;
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '#' || *p == '\n' || *p == '\0') continue;

    char kind[32];
    if (sscanf(p, "%31s", kind) != 1) continue;

    if (!strcmp(kind, "collection")) {
      char id[64], mode[16];
      if (sscanf(p, "%*s %63s %15s", id, mode) == 2)
        cur = pc_new_collection(s, id, !strcmp(mode, "absolute"));
    } else if (!strcmp(kind, "point")) {
      float z, y, x, w;
      if (cur && sscanf(p, "%*s %f %f %f %f", &z, &y, &x, &w) == 4)
        pc_add_point(cur, z, y, x, w);
    } else if (!strcmp(kind, "cannotlink") || !strcmp(kind, "mustlink")) {
      link_constraint lk;
      lk.cannot = !strcmp(kind, "cannotlink");
      if (sscanf(p, "%*s %f %f %f %f %f %f", &lk.az, &lk.ay, &lk.ax,
                 &lk.bz, &lk.by, &lk.bx) == 6)
        pc_add_link(s, lk);
    }
  }
  fclose(f);
  return 0;
}

void pc_free(pc_set *s) {
  for (int i = 0; i < s->ncols; i++) free(s->cols[i].pts);
  free(s->cols);
  free(s->links);
  memset(s, 0, sizeof *s);
}
