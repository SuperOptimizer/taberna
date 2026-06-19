/* test_features.c — cubical_features must return, per dim-1 tunnel, a valid
 * representative 1-cycle (every vertex incident to an even number of loop edges).
 * Feed the INVERTED mask (FG->0) so foreground tunnels are off-diagonal. */
#include <stdio.h>
#include <stdlib.h>
#include "topo/cubical.h"

int main(void){
  int nz=14,ny=14,nx=14; size_t N=(size_t)nz*ny*nx;
  f32 *f=malloc(N*sizeof(f32));
  // solid bar spanning z, with a through-tunnel drilled along z (b1=1)
  for(size_t i=0;i<N;i++) f[i]=1.0f;             // background = 1
  for(int z=0;z<nz;z++)for(int y=3;y<=10;y++)for(int x=3;x<=10;x++) f[((size_t)z*ny+y)*nx+x]=0.0f; // FG=0
  for(int z=0;z<nz;z++) f[((size_t)z*ny+6)*nx+6]=1.0f;  // drill tunnel
  int nf; pers_feat *ft=cubical_features(f,nz,ny,nx,&nf);
  int d1=0, badcycle=0;
  for(int i=0;i<nf;i++){
    if(ft[i].dim!=1 || ft[i].death>=TOPO_INF) continue;
    d1++;
    // verify even vertex degree: each edge gid = N + axis*N + v connects v and v+e_axis
    int *deg=calloc(N,sizeof(int));
    for(int t=0;t<ft[i].ncells;t++){
      long gid=ft[i].cells[t]; long v=gid%(long)N; int axis=(gid-(long)N)/(long)N;
      long v2 = v + (axis==0?1: axis==1?nx : nx*ny);
      deg[v]^=1; if(v2>=0&&v2<(long)N) deg[v2]^=1;
    }
    for(size_t v=0;v<N;v++) if(deg[v]&1){ badcycle=1; break; }
    free(deg);
  }
  printf("dim-1 off-diagonal features: %d (want 1);  cycle valid: %s\n",
         d1, badcycle?"NO":"yes");
  for(int i=0;i<nf;i++) free(ft[i].cells);
  free(ft); free(f);
  int ok = (d1==1 && !badcycle);
  printf(ok?"PASS\n":"FAIL\n");
  return ok?0:1;
}
