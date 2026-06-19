/* test_persist0.c — validate fast union-find dim-0 persistence against the generic
 * boundary-matrix engine, and sanity-check simplify_dim0. Self-contained (random
 * fields, fixed seed). Exit 0 = pass. */
#include <stdio.h>
#include <stdlib.h>
#include "topo/cubical.h"
#include "topo/persist0.h"

static int cmpf(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return (x>y)-(x<y); }

/* collect sorted dim-0 off-diagonal persistences from a pair array */
static int dim0_persist(const pers_pair*p,int n,float*out){
  int m=0; for(int i=0;i<n;i++) if(p[i].dim==0 && p[i].death<TOPO_INF) out[m++]=p[i].death-p[i].birth;
  qsort(out,m,sizeof(float),cmpf); return m;
}

int main(void){
  unsigned seed=12345; int fails=0;
  for(int trial=0; trial<200; trial++){
    int nz=3+seed%5; seed=seed*1103515245u+12345u;
    int ny=3+seed%5; seed=seed*1103515245u+12345u;
    int nx=3+seed%5; seed=seed*1103515245u+12345u;
    size_t N=(size_t)nz*ny*nx;
    f32*f=malloc(N*sizeof(f32));
    for(size_t i=0;i<N;i++){ seed=seed*1103515245u+12345u; f[i]=(f32)(seed%7); }
    int n1,n2; pers_pair*a=cubical_persistence(f,nz,ny,nx,&n1);
    pers_pair*b=dim0_persistence(f,nz,ny,nx,&n2);
    float*pa=malloc(N*sizeof(float)),*pb=malloc(N*sizeof(float));
    int ma=dim0_persist(a,n1,pa), mb=dim0_persist(b,n2,pb);
    int ok=(ma==mb); for(int i=0;i<ma&&ok;i++) if(pa[i]!=pb[i]) ok=0;
    if(!ok){ fails++; if(fails<=5) printf("FAIL %dx%dx%d: cubical=%d uf=%d dim0 pairs\n",nz,ny,nx,ma,mb); }
    free(a);free(b);free(pa);free(pb);free(f);
  }
  // simplify sanity: deep global min at idx0, shallow secondary min (7) that
  // merges into it at saddle 9 (persistence 2). tau=3 should raise it to 9; the
  // deeper min must survive.
  { int nz=1,ny=1,nx=4; f32 g[4]={0,9,7,9};
    simplify_dim0(g,nz,ny,nx,3.0f);
    if(g[2]!=9.0f){ printf("FAIL simplify: shallow basin not filled (%.1f)\n",g[2]); fails++; }
    if(g[0]!=0.0f){ printf("FAIL simplify: deep min wrongly raised (%.1f)\n",g[0]); fails++; }
  }
  // a basin DEEPER than tau must be preserved
  { int nz=1,ny=1,nx=4; f32 g[4]={0,9,2,9};
    simplify_dim0(g,nz,ny,nx,3.0f);   // persistence 7 > 3 => keep
    if(g[2]!=2.0f){ printf("FAIL simplify: deep basin wrongly filled (%.1f)\n",g[2]); fails++; }
  }
  printf(fails?"\n%d FAILED\n":"\nALL PASS (dim0 union-find == generic engine, simplify ok)\n",fails);
  return fails?1:0;
}
