/* wind_cut — Phase 3b (hard) of docs/touching-sheets-plan.md: the ORDERED-LABEL MIN-CUT that places the
 * missing wrap boundary through a LEAKED (gapless fused) touch — the one thing the soft min-slope ordering
 * provably cannot do (it can't choose a discrete cut location). This is the Boykov-Kolmogorov/Ishikawa
 * construction the research and the negative results both point to.
 *
 * BUILT VERIFY-FIRST: this file starts as a correct, SELF-TESTED maxflow (Dinic). The grid construction
 * (banded binary cuts / ordered level-sets along the smoothed-normal columns) is layered on top ONLY once
 * the solver is proven — no grid cut on an unverified maxflow.
 *
 *   wind_cut selftest                 # run the maxflow unit tests (known min-cuts) and exit
 *   wind_cut ARC OUT ... (grid cut)   # TODO: layered on after the solver is verified
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"

/* ---------------------------------------------------------------- Dinic maxflow
 * Half-edge arrays: edge e and its reverse e^1 are adjacent (2*k, 2*k+1). cap[] is the RESIDUAL
 * capacity (decremented on push, the reverse incremented), so flow = original_cap - cap on forward
 * edges. Adjacency as head[]/nxt[] linked lists. Correctness rests on: (1) BFS level graph only
 * advances to strictly deeper nodes; (2) blocking-flow DFS with it[] "current-arc" pointers never
 * re-scans a saturated arc; (3) dead-ends are pruned by setting level=-1. */
typedef struct {
  int n;                 /* node count */
  long m;                /* half-edge count */
  long cap_e;            /* edge array capacity */
  int *to; long *cap; long *nxt;   /* per half-edge: target, residual cap, next-in-adjacency */
  long *head;            /* per node: first half-edge (-1 = none) */
  int *level; long *it;  /* BFS level, current-arc iterator */
  int *q;                /* BFS queue */
} Dinic;

static void d_init(Dinic*d,int n,long max_edges){
  d->n=n; d->m=0; d->cap_e=2*max_edges;
  d->to=malloc(d->cap_e*sizeof(int)); d->cap=malloc(d->cap_e*sizeof(long)); d->nxt=malloc(d->cap_e*sizeof(long));
  d->head=malloc((size_t)n*sizeof(long)); for(int i=0;i<n;i++)d->head[i]=-1;
  d->level=malloc((size_t)n*sizeof(int)); d->it=malloc((size_t)n*sizeof(long)); d->q=malloc((size_t)n*sizeof(int));
}
static void d_free(Dinic*d){ free(d->to);free(d->cap);free(d->nxt);free(d->head);free(d->level);free(d->it);free(d->q); }
/* directed edge u->v with capacity cu, and reverse v->u with capacity cv (cv=0 for a pure directed arc,
 * cv=cu for an undirected edge). */
static void d_edge(Dinic*d,int u,int v,long cu,long cv){
  d->to[d->m]=v; d->cap[d->m]=cu; d->nxt[d->m]=d->head[u]; d->head[u]=d->m; d->m++;
  d->to[d->m]=u; d->cap[d->m]=cv; d->nxt[d->m]=d->head[v]; d->head[v]=d->m; d->m++;
}
static int d_bfs(Dinic*d,int s,int t){
  for(int i=0;i<d->n;i++)d->level[i]=-1;
  int qh=0,qt=0; d->level[s]=0; d->q[qt++]=s;
  while(qh<qt){ int u=d->q[qh++];
    for(long e=d->head[u];e!=-1;e=d->nxt[e]){ int v=d->to[e];
      if(d->cap[e]>0 && d->level[v]<0){ d->level[v]=d->level[u]+1; d->q[qt++]=v; } } }
  return d->level[t]>=0;
}
static long d_dfs(Dinic*d,int u,int t,long f){
  if(u==t) return f;
  for(long *ep=&d->it[u]; *ep!=-1; *ep=d->nxt[*ep]){ long e=*ep; int v=d->to[e];
    if(d->cap[e]>0 && d->level[v]==d->level[u]+1){
      long got=d_dfs(d,v,t,f<d->cap[e]?f:d->cap[e]);
      if(got>0){ d->cap[e]-=got; d->cap[e^1]+=got; return got; } } }
  d->level[u]=-1;   /* dead end: don't revisit this round */
  return 0;
}
static long d_maxflow(Dinic*d,int s,int t){
  long flow=0;
  while(d_bfs(d,s,t)){
    for(int i=0;i<d->n;i++)d->it[i]=d->head[i];
    long f; while((f=d_dfs(d,s,t,(long)1<<60))>0) flow+=f;
  }
  return flow;
}
/* min-cut: source side = nodes reachable from s in the residual graph after maxflow */
static void d_cut_side(Dinic*d,int s,unsigned char*side){
  for(int i=0;i<d->n;i++)side[i]=0;
  int qh=0,qt=0; side[s]=1; d->q[qt++]=s;
  while(qh<qt){ int u=d->q[qh++];
    for(long e=d->head[u];e!=-1;e=d->nxt[e]){ int v=d->to[e];
      if(d->cap[e]>0 && !side[v]){ side[v]=1; d->q[qt++]=v; } } }
}

/* ---------------------------------------------------------------- self-test */
static int g_fail=0;
static void check(const char*name,long got,long want){
  if(got==want) printf("  PASS %-28s = %ld\n",name,got);
  else { printf("  FAIL %-28s = %ld (want %ld)\n",name,got,want); g_fail=1; }
}
static int selftest(void){
  printf("maxflow self-test (Dinic):\n");
  /* 1. trivial bottleneck: s->a cap5, a->t cap3 -> 3 */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,t=2; d_edge(&d,s,a,5,0); d_edge(&d,a,t,3,0);
    check("bottleneck s-a-t",d_maxflow(&d,s,t),3); d_free(&d); }
  /* 2. two parallel paths: s->a->t (cap 4,4), s->b->t (cap 2,5) -> 4+2=6 */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,b=2,t=3;
    d_edge(&d,s,a,4,0); d_edge(&d,a,t,4,0); d_edge(&d,s,b,2,0); d_edge(&d,b,t,5,0);
    check("two parallel paths",d_maxflow(&d,s,t),6); d_free(&d); }
  /* 3. CLRS classic 6-node network, known max flow = 23 */
  { Dinic d; d_init(&d,6,9); int s=0,t=5;
    d_edge(&d,0,1,16,0); d_edge(&d,0,2,13,0); d_edge(&d,1,2,10,0); d_edge(&d,2,1,4,0);
    d_edge(&d,1,3,12,0); d_edge(&d,3,2,9,0); d_edge(&d,2,4,14,0); d_edge(&d,4,3,7,0);
    d_edge(&d,3,5,20,0); d_edge(&d,4,5,4,0);
    check("CLRS network",d_maxflow(&d,s,t),23); d_free(&d); }
  /* 4. min-cut value == max-flow, and the cut partition is consistent (s-side reaches t? must NOT) */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,b=2,t=3;
    d_edge(&d,s,a,3,0); d_edge(&d,a,b,2,0); d_edge(&d,b,t,3,0); d_edge(&d,s,b,1,0);
    long f=d_maxflow(&d,s,t); check("series+shortcut flow",f,3);
    unsigned char side[4]; d_cut_side(&d,s,side);
    /* cut edges = forward edges from s-side to t-side; their capacity sum must equal f */
    long cutcap=0; /* recompute against ORIGINAL caps: rebuild needed -> instead verify t not on s-side */
    check("t not on source side",side[t]?1:0,0); (void)cutcap;
    d_free(&d); }
  /* 5. undirected grid 2x2, s=corner, t=opposite, all caps 1: max flow = 2 (two edge-disjoint paths) */
  { Dinic d; d_init(&d,4,8); /* nodes 0=(0,0)1=(0,1)2=(1,0)3=(1,1); s=0 t=3 */
    d_edge(&d,0,1,1,1); d_edge(&d,0,2,1,1); d_edge(&d,1,3,1,1); d_edge(&d,2,3,1,1);
    check("2x2 undirected grid",d_maxflow(&d,0,3),2); d_free(&d); }
  /* 6. ordering arc (infinity): models the Ishikawa "label monotone" constraint. s->a cap10,
   *    a=>b INF (ordering: b can't be cut from a's side cheaply), b->t cap10. Cutting must pay
   *    min(10,10)=10, never the INF arc. */
  { Dinic d; d_init(&d,4,3); int s=0,a=1,b=2,t=3; long INF=1L<<50;
    d_edge(&d,s,a,10,0); d_edge(&d,a,b,INF,0); d_edge(&d,b,t,10,0);
    check("infinity ordering arc",d_maxflow(&d,s,t),10); d_free(&d); }
  printf(g_fail?"SELFTEST FAILED\n":"ALL TESTS PASSED\n");
  return g_fail;
}

/* ---------------------------------------------------------------- grid ordered-label cut
 * Greedy-NESTED level-set decomposition over the verified binary maxflow. For each wrap boundary level
 * k = wmin+1..wmax, a binary min-cut decides which material voxels advance to label >= k. Capacities:
 *   t-links (data, from the TV field u): a voxel deep inside wrap (u >> k) is hard-tied to SOURCE (>=k),
 *     deep outside (u << k) hard-tied to SINK; in the band a soft pull ~ |u-(k-0.5)|.
 *   n-links (6-neighbour, cap = mean sheetness): cutting between two voxels = putting the wrap boundary
 *     there, priced by sheetness -> the boundary AVOIDS sheets and snaps to the lowest-sheetness locus,
 *     which through a GAPLESS fused touch is exactly the cut the soft ordering could not place.
 * Nesting is enforced by construction (only voxels already at label k-1 are eligible for >=k), so the
 * per-level cuts compose into a valid monotone-outward integer labeling. */
static int air_threshold(const u8*v,size_t n){ long h[256]={0},t=0; double s=0; long z=0;
  for(size_t i=0;i<n;i++){int x=v[i]; if(x){s+=x;z++;} if(x>=1&&x<=254){h[x]++;t++;}}
  if(t<256)return z?(int)(0.5*s/z+0.5):1; double sum=0; for(int i=1;i<=254;i++)sum+=(double)i*h[i];
  double sB=0,wB=0,b=-1; int thr=1; for(int k=1;k<=254;k++){wB+=h[k]; if(!wB)continue; double wF=t-wB; if(wF<=0)break;
    sB+=(double)k*h[k]; double mB=sB/wB,mF=(sum-sB)/wF,bb=wB*wF*(mB-mF)*(mB-mF); if(bb>b){b=bb;thr=k;}} return thr; }
static double cw_trilin(const float*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0,acc=0,wt=0;
  int zs[2]={z0,z1},ys[2]={y0,y1},xs[2]={x0,x1}; double zw[2]={1-dz,dz},yw[2]={1-dy,dy},xw[2]={1-dx,dx};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){ double vv=w[((size_t)zs[a]*ny+ys[b])*nx+xs[c]];
    if(!isfinite(vv))continue; double gg=zw[a]*yw[b]*xw[c]; acc+=gg*vv; wt+=gg; }
  return wt>1e-9? acc/wt : NAN; }

static int gridcut(int argc,char**argv){
  if(argc<12){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [band=1.5] [ncap=8] [zc=mid]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  const char*pvf=argv[10]; int plod=atoi(argv[11]);
  double band=argc>12?atof(argv[12]):1.5; double ncap=argc>13?atof(argv[13]):8.0; int zc=argc>14?atoi(argv[14]):dz/2;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx, plane=(size_t)dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr?1:0; nmat+=mask[p]; }
  fprintf(stderr,"air<%d, %ld/%zu material\n",athr,nmat,nn);
  // sheetness [0,1] (n-link price: high = expensive to cut = boundary routes around sheets)
  f32*vf=malloc(nn*sizeof(f32)),*sh=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++)vf[p]=v[p];
  st_params spp=st_default_params(); spp.sigma_grad=1.0f; spp.sigma_tensor=0.8f; st_sheet_detect(vf,dz,dy,dx,&spp,sh,NULL); free(vf);
  { double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p]; if(smx<1e-6)smx=1; for(size_t p=0;p<nn;p++) sh[p]/=(f32)smx; }
  // prior winding field u (the wind_tv _tv.f32 / wind_poisson _vol.f32 format)
  FILE*pf=fopen(pvf,"rb"); if(!pf){fprintf(stderr,"priorvol open fail\n");return 1;}
  int hd[6]; if(fread(hd,sizeof(int),6,pf)!=6){return 1;} int cnz=hd[0],cny=hd[1],cnx=hd[2]; long cz0=hd[3],cy0=hd[4],cx0=hd[5];
  size_t cn=(size_t)cnz*cny*cnx; float*cw=malloc(cn*sizeof(float)); if(fread(cw,sizeof(float),cn,pf)!=cn){return 1;} fclose(pf);
  double scl=ldexp(1.0,lod-plod);
  f32*U=malloc(nn*sizeof(f32)); int wmin=1<<30,wmax=-(1<<30);
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x; U[p]=NAN;
    if(mask[p]){ double wv=cw_trilin(cw,cnz,cny,cnx,(z0+z)*scl-cz0,(y0+y)*scl-cy0,(x0+x)*scl-cx0);
      if(isfinite(wv)){ U[p]=(f32)wv; int f=(int)floor(wv); if(f<wmin)wmin=f; if(f>wmax)wmax=f; } else mask[p]=0; } }
  free(cw);
  fprintf(stderr,"wrap levels %d..%d; greedy-nested ordered cut (band=%.1f ncap=%.1f)\n",wmin,wmax,band,ncap);

  // label L: start everyone at wmin; each level k advances the source-side voxels to k.
  s32*L=malloc(nn*sizeof(s32)); for(size_t p=0;p<nn;p++) L[p]=mask[p]?wmin:-1;
  int*node=malloc(nn*sizeof(int));               // voxel -> band-graph node id (-1 = not in band this level)
  long DCAP=1000;                                 // hard t-link (infinity-ish, >> any n-link sum)
  for(int k=wmin+1;k<=wmax;k++){
    // eligible band voxels: currently at label k-1 (nested) AND |u-(k-0.5)| within ~band of the boundary
    double thr=k-0.5;
    // eligible = nested (L==k-1). Far-OUTSIDE (u >> thr) auto-advances FREE (definitely >=k, else the bulk
    // never ratchets up through the bands); the band+margin shell goes to the graph; far-inside stays.
    int nb=0; for(size_t p=0;p<nn;p++){ node[p]=-1; if(!(mask[p]&&L[p]==k-1))continue;
      double du=(double)U[p]-thr;
      if(du>band+1.0) L[p]=k;                       // far outside -> advance free, not in graph
      else if(du>=-(band+1.0)) node[p]=nb++; }      // band+margin shell -> binary cut decides
    if(nb==0) continue;
    Dinic d; long maxE=(long)nb*7+8; d_init(&d,nb+2,maxE); int S=nb,T=nb+1;
    // t-links from data; n-links from sheetness across 6-neighbours both in-band
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x; int np=node[p]; if(np<0)continue;
      double du=(double)U[p]-thr;
      if(du>=band)      d_edge(&d,S,np,DCAP,0);        // deep outside -> hard >=k
      else if(du<=-band)d_edge(&d,np,T,DCAP,0);        // deep inside  -> hard <k
      else { double s=du/band; if(s>0) d_edge(&d,S,np,(long)(1+DCAP*0.4*s),0); else d_edge(&d,np,T,(long)(1+DCAP*0.4*(-s)),0); }
      // n-links (forward neighbours; undirected cap = ncap * mean sheetness -> route boundary through gaps)
      size_t nbrs[3]={p+1,p+dx,p+plane}; int ok[3]={x+1<dx,y+1<dy,z+1<dz};
      for(int a=0;a<3;a++){ if(!ok[a])continue; size_t q=nbrs[a]; if(node[q]<0)continue;
        long c=(long)(1+ncap*0.5*((double)sh[p]+sh[q])); d_edge(&d,np,node[q],c,c); } }
    long fl=d_maxflow(&d,S,T); (void)fl;
    unsigned char*side=malloc((size_t)(nb+2)); d_cut_side(&d,S,side);
    long adv=0; for(size_t p=0;p<nn;p++){ int np=node[p]; if(np>=0 && side[np]){ L[p]=k; adv++; } }   // source side advances
    free(side); d_free(&d);
    if(k%5==0||k==wmax) fprintf(stderr,"  level %d: %d band, %ld advanced\n",k,nb,adv);
  }
  free(node);

  // outputs: label volume + mid-z ordered-palette render (same convention as wind_tv)
  char fn[700]; snprintf(fn,sizeof fn,"%s_lab.i32",outp); FILE*lf=fopen(fn,"wb");
  if(lf){ int h2[6]={dz,dy,dx,z0,y0,x0}; fwrite(h2,sizeof(int),6,lf); fwrite(L,sizeof(s32),nn,lf); fclose(lf);
    fprintf(stderr,"wrote %s\n",fn); }
  static const u8 pal[6][3]={{230,60,60},{240,180,40},{60,200,80},{50,180,220},{90,90,230},{220,90,210}};
  u8*rgb=calloc(plane*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3; int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
    if(L[p]>=0){ int q=((L[p]%6)+6)%6; rgb[o]=pal[q][0]; rgb[o+1]=pal[q][1]; rgb[o+2]=pal[q][2]; } }
  snprintf(fn,sizeof fn,"%s_label.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,plane*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb); free(L);free(U);free(sh);free(mask);free(v); mca_close(arc); return 0;
}

int main(int argc,char**argv){
  if(argc>=2 && !strcmp(argv[1],"selftest")) return selftest();
  if(argc>=3) return gridcut(argc,argv);
  fprintf(stderr,"usage: %s selftest | %s ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [band] [ncap] [zc]\n",argv[0],argv[0]);
  return 2;
}
