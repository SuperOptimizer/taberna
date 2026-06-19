/* topo.c — see topo.h. Exact discrete topology via O(N) cell counting + labeling. */
#include "eval/topo.h"

#include <stdlib.h>
#include <string.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

/* ---- union-find ----------------------------------------------------------- */
typedef struct { u32 *parent; } uf;
static u32 uf_find(u32 *p, u32 a) {
  while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; }
  return a;
}
static void uf_union(u32 *p, u32 a, u32 b) {
  a = uf_find(p, a); b = uf_find(p, b);
  if (a != b) p[a < b ? b : a] = (a < b ? a : b);  // point larger root at smaller
}

/* Connected components by union-find over already-seen lower neighbors. conn=6
 * uses the 3 face neighbors (-z,-y,-x); conn=26 adds the 9+3 lower diagonal ones.
 * Returns component count; fills `labels` with compacted 1-based ids. */
u32 cc_label(const u8 *mask, int nz, int ny, int nx, topo_conn conn, u32 *labels) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  // provisional per-voxel id = its own index+1 for fg, 0 for bg
  u32 *p = (u32 *)malloc((n + 1) * sizeof(u32));
  for (size_t i = 0; i <= n; i++) p[i] = (u32)i;

  // lower-neighbor offsets (must have strictly smaller scan order)
  static const int off6[3][3]  = {{-1,0,0},{0,-1,0},{0,0,-1}};
  static const int off26[13][3] = {
    {-1,0,0},{0,-1,0},{0,0,-1},
    {-1,-1,0},{-1,1,0},{-1,0,-1},{-1,0,1},{0,-1,-1},{0,-1,1},
    {-1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1}};
  const int (*off)[3] = conn == TOPO_CONN26 ? off26 : off6;
  int noff = conn == TOPO_CONN26 ? 13 : 3;

  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (!mask[i]) continue;
        for (int k = 0; k < noff; k++) {
          int zz = z + off[k][0], yy = y + off[k][1], xx = x + off[k][2];
          if (zz < 0 || yy < 0 || yy >= ny || xx < 0 || xx >= nx) continue;
          size_t j = IDX(zz, yy, xx);
          if (mask[j]) uf_union(p, (u32)(i + 1), (u32)(j + 1));
        }
      }

  // compact roots -> 1-based labels
  u32 *remap = (u32 *)calloc(n + 1, sizeof(u32));
  u32 count = 0;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (!mask[i]) { labels[i] = 0; continue; }
        u32 r = uf_find(p, (u32)(i + 1));
        if (!remap[r]) remap[r] = ++count;
        labels[i] = remap[r];
      }
  free(p); free(remap);
  return count;
}

/* ---- exact Euler characteristic of the foreground cubical complex ---------
 * Treat each fg voxel as a closed unit cube. chi = V - E + F - C over the DISTINCT
 * cells of the union. A k-cell exists iff at least one of the voxels incident to it
 * is foreground; counting each cell once over the (nz+1)x(ny+1)x(nx+1) lattice of
 * vertices / the three edge lattices / the three face lattices gives exact V,E,F. */
long euler_characteristic(const u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  long V = 0, E = 0, F = 0, C = 0;

  // helper: is voxel fg (with bounds check returning 0 outside)
  #define FG(z, y, x) (((z) >= 0 && (z) < nz && (y) >= 0 && (y) < ny && \
                        (x) >= 0 && (x) < nx) ? (mask[IDX(z, y, x)] != 0) : 0)

  // Cubes: one per fg voxel.
  for (size_t i = 0, n = (size_t)nz * nynx; i < n; i++) C += (mask[i] != 0);

  // Vertices: lattice point (z,y,x) in [0..nz]x[0..ny]x[0..nx]; incident to the 8
  // voxels (z-1..z, y-1..y, x-1..x). Exists if any is fg.
  for (int z = 0; z <= nz; z++)
    for (int y = 0; y <= ny; y++)
      for (int x = 0; x <= nx; x++) {
        int any = FG(z-1,y-1,x-1)||FG(z-1,y-1,x)||FG(z-1,y,x-1)||FG(z-1,y,x)||
                  FG(z,y-1,x-1)||FG(z,y-1,x)||FG(z,y,x-1)||FG(z,y,x);
        V += any;
      }

  // Edges: 3 orientations. An edge along axis A at the appropriate lattice is
  // incident to the 4 voxels sharing it. Exists if any of the 4 is fg.
  // x-edges: span x->x+1 at vertex (z,y,*); incident voxels (z-1..z, y-1..y, x).
  for (int z = 0; z <= nz; z++)
    for (int y = 0; y <= ny; y++)
      for (int x = 0; x < nx; x++)
        E += (FG(z-1,y-1,x)||FG(z-1,y,x)||FG(z,y-1,x)||FG(z,y,x)) != 0;
  // y-edges: incident voxels (z-1..z, y, x-1..x).
  for (int z = 0; z <= nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x <= nx; x++)
        E += (FG(z-1,y,x-1)||FG(z-1,y,x)||FG(z,y,x-1)||FG(z,y,x)) != 0;
  // z-edges: incident voxels (z, y-1..y, x-1..x).
  for (int z = 0; z < nz; z++)
    for (int y = 0; y <= ny; y++)
      for (int x = 0; x <= nx; x++)
        E += (FG(z,y-1,x-1)||FG(z,y-1,x)||FG(z,y,x-1)||FG(z,y,x)) != 0;

  // Faces: 3 orientations, each shared by 2 voxels across it.
  // xy-faces (normal z) at z-plane: voxels (z-1, y, x) and (z, y, x).
  for (int z = 0; z <= nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++)
        F += (FG(z-1,y,x)||FG(z,y,x)) != 0;
  // xz-faces (normal y): voxels (z, y-1, x) and (z, y, x).
  for (int z = 0; z < nz; z++)
    for (int y = 0; y <= ny; y++)
      for (int x = 0; x < nx; x++)
        F += (FG(z,y-1,x)||FG(z,y,x)) != 0;
  // yz-faces (normal x): voxels (z, y, x-1) and (z, y, x).
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x <= nx; x++)
        F += (FG(z,y,x-1)||FG(z,y,x)) != 0;

  #undef FG
  return V - E + F - C;
}

/* ---- Betti numbers --------------------------------------------------------
 * b0 = #(26-conn fg components). b2 = #(6-conn bg components fully enclosed, i.e.
 * not touching the volume border). b1 = b0 + b2 - chi. */
topo_betti betti_numbers(const u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  topo_betti r = {0, 0, 0, 0};

  u32 *lab = (u32 *)malloc(n * sizeof(u32));
  r.b0 = (long)cc_label(mask, nz, ny, nx, TOPO_CONN26, lab);

  // background mask, label its 6-conn components, count those not touching border.
  u8 *bg = (u8 *)malloc(n);
  for (size_t i = 0; i < n; i++) bg[i] = mask[i] ? 0 : 1;
  u32 nbg = cc_label(bg, nz, ny, nx, TOPO_CONN6, lab);
  if (nbg) {
    u8 *touches = (u8 *)calloc(nbg + 1, 1);
    for (int z = 0; z < nz; z++)
      for (int y = 0; y < ny; y++)
        for (int x = 0; x < nx; x++) {
          if (z && z != nz-1 && y && y != ny-1 && x && x != nx-1) continue; // border only
          u32 l = lab[IDX(z, y, x)];
          if (l) touches[l] = 1;
        }
    long enclosed = 0;
    for (u32 c = 1; c <= nbg; c++) if (!touches[c]) enclosed++;
    r.b2 = enclosed;
    free(touches);
  }
  free(bg); free(lab);

  r.chi = euler_characteristic(mask, nz, ny, nx);
  r.b1 = r.b0 + r.b2 - r.chi;
  return r;
}

/* Euler characteristic of the (6-conn fg) voxel-as-vertex cubical complex:
 * chi = V - E + F - C  (vertices, 6-adjacency edges, unit squares, unit cubes). */
static long euler_characteristic_6(const u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  long V = 0, E = 0, F = 0, C = 0;
  #define M(z,y,x) (mask[IDX(z,y,x)] != 0)
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        if (!M(z,y,x)) continue;
        V++;
        if (x+1<nx && M(z,y,x+1)) E++;
        if (y+1<ny && M(z,y+1,x)) E++;
        if (z+1<nz && M(z+1,y,x)) E++;
        if (x+1<nx && y+1<ny && M(z,y,x+1) && M(z,y+1,x) && M(z,y+1,x+1)) F++;          // xy
        if (x+1<nx && z+1<nz && M(z,y,x+1) && M(z+1,y,x) && M(z+1,y,x+1)) F++;          // xz
        if (y+1<ny && z+1<nz && M(z,y+1,x) && M(z+1,y,x) && M(z+1,y+1,x)) F++;          // yz
        if (x+1<nx && y+1<ny && z+1<nz &&
            M(z,y,x+1)&&M(z,y+1,x)&&M(z,y+1,x+1)&&M(z+1,y,x)&&M(z+1,y,x+1)&&M(z+1,y+1,x)&&M(z+1,y+1,x+1)) C++;
      }
  #undef M
  return V - E + F - C;
}

topo_betti betti_numbers_6(const u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  topo_betti r = {0,0,0,0};
  u32 *lab = (u32 *)malloc(n * sizeof(u32));
  r.b0 = (long)cc_label(mask, nz, ny, nx, TOPO_CONN6, lab);    // 6-conn foreground
  u8 *bg = (u8 *)malloc(n);
  for (size_t i = 0; i < n; i++) bg[i] = mask[i] ? 0 : 1;
  u32 nbg = cc_label(bg, nz, ny, nx, TOPO_CONN26, lab);        // 26-conn background
  if (nbg) {
    u8 *touch = (u8 *)calloc(nbg + 1, 1);
    for (int z = 0; z < nz; z++)
      for (int y = 0; y < ny; y++)
        for (int x = 0; x < nx; x++) {
          if (z && z!=nz-1 && y && y!=ny-1 && x && x!=nx-1) continue;
          u32 l = lab[IDX(z,y,x)]; if (l) touch[l] = 1;
        }
    for (u32 c = 1; c <= nbg; c++) if (!touch[c]) r.b2++;
    free(touch);
  }
  free(bg); free(lab);
  r.chi = euler_characteristic_6(mask, nz, ny, nx);
  r.b1 = r.b0 + r.b2 - r.chi;
  return r;
}

long region_b1(const u8 *mask, int nz, int ny, int nx,
               int z0, int y0, int x0, int dz, int dy, int dx) {
  if (z0 < 0) { dz += z0; z0 = 0; }
  if (y0 < 0) { dy += y0; y0 = 0; }
  if (x0 < 0) { dx += x0; x0 = 0; }
  if (z0 + dz > nz) dz = nz - z0;
  if (y0 + dy > ny) dy = ny - y0;
  if (x0 + dx > nx) dx = nx - x0;
  if (dz <= 0 || dy <= 0 || dx <= 0) return 0;

  size_t nynx = (size_t)ny * nx;
  u8 *sub = (u8 *)malloc((size_t)dz * dy * dx);
  for (int z = 0; z < dz; z++)
    for (int y = 0; y < dy; y++)
      for (int x = 0; x < dx; x++)
        sub[((size_t)z * dy + y) * dx + x] = mask[IDX(z0 + z, y0 + y, x0 + x)];
  topo_betti b = betti_numbers(sub, dz, dy, dx);
  free(sub);
  return b.b1;
}
