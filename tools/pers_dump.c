/* pers_dump.c — dump the persistence diagram of a small raw f32 volume, for
 * validating taberna's cubical PH against the Betti-Matching-3D oracle.
 * Reads a headerless little-endian f32 volume (nz,ny,nx from argv). Prints
 * off-diagonal pairs "dim birth death" (one per line), sorted, then essentials.
 *
 * Usage: pers_dump FIELD.f32 nz ny nx
 */
#include <stdio.h>
#include <stdlib.h>
#include "topo/cubical.h"

static int cmp(const void*a,const void*b){
  const pers_pair*x=a,*y=b;
  if(x->dim!=y->dim) return x->dim-y->dim;
  if(x->birth<y->birth) return -1; if(x->birth>y->birth) return 1;
  if(x->death<y->death) return -1; if(x->death>y->death) return 1; return 0;
}

int main(int argc,char**argv){
  if(argc!=5){fprintf(stderr,"usage: %s FIELD.f32 nz ny nx\n",argv[0]);return 2;}
  int nz=atoi(argv[2]),ny=atoi(argv[3]),nx=atoi(argv[4]);
  size_t n=(size_t)nz*ny*nx;
  f32*f=malloc(n*sizeof(f32));
  FILE*fp=fopen(argv[1],"rb"); if(!fp||fread(f,sizeof(f32),n,fp)!=n){fprintf(stderr,"read fail\n");return 1;} fclose(fp);
  int np; pers_pair*p=cubical_persistence(f,nz,ny,nx,&np);
  qsort(p,np,sizeof(pers_pair),cmp);
  for(int i=0;i<np;i++){
    if(p[i].death>=TOPO_INF) printf("ess dim %d birth %.4f\n",p[i].dim,p[i].birth);
    else printf("dim %d  %.4f -> %.4f\n",p[i].dim,p[i].birth,p[i].death);
  }
  free(p);free(f);
  return 0;
}
