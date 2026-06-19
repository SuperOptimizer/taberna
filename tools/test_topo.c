/* test_topo.c — self-contained unit tests for the topology + postproc toolkit.
 * Verifies exact discrete invariants on synthetic shapes (ball/cavity/tunnel) and
 * that the classical cleanup ops do what the Kaggle writeups claim:
 *   - majority (median) filter removes a thin bridge => unmerges components
 *   - remove_small_components drops dust
 *   - fill_holes kills a cavity (b2 -> 0)
 *   - plug_pinholes fills a 1-voxel hole
 * Exit code 0 = all pass. Run: tools/test_topo
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eval/topo.h"
#include "postproc/morph.h"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { \
  printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); g_fail++; } \
  else { printf("ok:   "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static u8 *vol(int nz, int ny, int nx) { return (u8 *)calloc((size_t)nz*ny*nx, 1); }
#define AT(m, z, y, x) (m)[((size_t)(z)*ny + (y))*nx + (x)]

static void fill_box(u8 *m, int ny, int nx, int z0,int z1,int y0,int y1,int x0,int x1) {
  for (int z=z0; z<=z1; z++) for (int y=y0; y<=y1; y++) for (int x=x0; x<=x1; x++)
    AT(m,z,y,x) = 1;
}

int main(void) {
  // --- 1. solid block: chi=1, b0=1, b1=0, b2=0 -------------------------------
  {
    int nz=12,ny=12,nx=12; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 2,9, 2,9, 2,9);
    topo_betti b = betti_numbers(m,nz,ny,nx);
    CHECK(b.b0==1 && b.b1==0 && b.b2==0 && b.chi==1,
          "solid block: b0=%ld b1=%ld b2=%ld chi=%ld (want 1,0,0,1)", b.b0,b.b1,b.b2,b.chi);
    free(m);
  }
  // --- 2. hollow block (1 cavity): chi=2, b2=1 -------------------------------
  {
    int nz=14,ny=14,nx=14; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 2,11, 2,11, 2,11);          // solid
    fill_box(m,ny,nx, 5,8, 5,8, 5,8);             // ensure solid then carve:
    for (int z=5; z<=8; z++) for (int y=5;y<=8;y++) for (int x=5;x<=8;x++) AT(m,z,y,x)=0;
    topo_betti b = betti_numbers(m,nz,ny,nx);
    CHECK(b.b0==1 && b.b1==0 && b.b2==1 && b.chi==2,
          "hollow block: b0=%ld b1=%ld b2=%ld chi=%ld (want 1,0,1,2)", b.b0,b.b1,b.b2,b.chi);
    // fill_holes should kill the cavity
    fill_holes(m,nz,ny,nx);
    topo_betti b2 = betti_numbers(m,nz,ny,nx);
    CHECK(b2.b2==0 && b2.b0==1, "fill_holes: b2=%ld b0=%ld (want 0,1)", b2.b2, b2.b0);
    free(m);
  }
  // --- 3. block with a through-tunnel (solid torus): chi=0, b1=1 -------------
  {
    int nz=14,ny=14,nx=14; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 0,nz-1, 3,10, 3,10);        // bar spanning all z
    for (int z=0; z<nz; z++) AT(m,z,6,6)=0;       // drill a z-through hole
    topo_betti b = betti_numbers(m,nz,ny,nx);
    CHECK(b.b0==1 && b.b1==1 && b.b2==0 && b.chi==0,
          "through-tunnel: b0=%ld b1=%ld b2=%ld chi=%ld (want 1,1,0,0)", b.b0,b.b1,b.b2,b.chi);
    free(m);
  }
  // --- 4. two blocks joined by a thin bridge; majority filter unmerges -------
  {
    int nz=20,ny=20,nx=40; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 5,14, 5,14, 2,11);          // block A
    fill_box(m,ny,nx, 5,14, 5,14, 28,37);         // block B
    for (int x=12; x<=27; x++) AT(m,10,10,x)=1;   // 1-voxel bridge
    topo_betti before = betti_numbers(m,nz,ny,nx);
    CHECK(before.b0==1, "bridged blocks before: b0=%ld (want 1)", before.b0);
    u8 *out=vol(nz,ny,nx);
    majority_filter(m,out,nz,ny,nx, 2);
    topo_betti after = betti_numbers(out,nz,ny,nx);
    CHECK(after.b0==2, "majority filter unmerge: b0=%ld (want 2)", after.b0);
    free(m); free(out);
  }
  // --- 5. dust removal -------------------------------------------------------
  {
    int nz=16,ny=16,nx=16; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 2,11, 2,11, 2,11);          // big component (~1000 vox)
    AT(m,14,14,14)=1; AT(m,14,14,13)=1;           // 2-voxel dust
    size_t rm = remove_small_components(m,nz,ny,nx, TOPO_CONN26, 100);
    topo_betti b = betti_numbers(m,nz,ny,nx);
    CHECK(rm==1 && b.b0==1, "dust removal: removed=%zu b0=%ld (want 1,1)", rm, b.b0);
    free(m);
  }
  // --- 6. pinhole plug -------------------------------------------------------
  {
    int nz=9,ny=9,nx=9; u8 *m=vol(nz,ny,nx);
    fill_box(m,ny,nx, 2,6, 2,6, 2,6);
    AT(m,4,4,4)=0;                                // single interior hole
    size_t p = plug_pinholes(m,nz,ny,nx);
    CHECK(p==1 && AT(m,4,4,4)==1, "plug_pinholes: plugged=%zu (want 1)", p);
    free(m);
  }

  printf(g_fail ? "\n%d CHECK(s) FAILED\n" : "\nALL TESTS PASSED\n", g_fail);
  return g_fail ? 1 : 0;
}
