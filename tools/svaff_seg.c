/* svaff_seg — SUPERVOXEL signed-affinity + Mutex-Watershed, the classical touching-sheets segmenter.
 *
 * The voxel-resolution path (affinity_seg) over-fragments: the signed repulsion (across-sheet contact
 * through high-sheetness material = crossing into the touching next wrap) is real, but at voxel scale
 * the normals are noisy and a single face's across*sheet can't be trusted -> 25k shards. The designed
 * fix (snic.h: "supervoxel coarsening is the path to TB scale, slots in behind the same sgraph"): run
 * SNIC supervoxels FIRST (they hug the sheets), then build the SAME signed affinity at SUPERVOXEL
 * resolution, AGGREGATING the across/sheetness signal over the whole contact patch between two
 * supervoxels (where the weak per-face signal averages into a decision) before MWS partitions it.
 *
 * Two pieces this tool adds over the existing library (the missing bits snic.h/affinity.h flagged):
 *   1. per-supervoxel ORIENTATION  : dominant eigenvector of the sheetness-weighted tensor sum
 *                                    Sum w*(n n^T) over the supervoxel's voxels (sign-ambiguity-safe,
 *                                    same definition as st_sheet_detect's normal = eigvec of largest l).
 *   2. supervoxel-level AFFINITY    : for each adjacent supervoxel pair, mean over all boundary faces of
 *                                    across=0.5(|Ni.d|+|Nj.d|) and boundary sheetness -> one signed edge.
 *
 *   svaff_seg ARC OUT lod z0 y0 x0 dz dy dx [d_seed=4] [compact=20] [krepel=1.0] [kattract=1.0] [zc=mid]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/snic.h"
#include "segmentation/affinity.h"
#include "segmentation/partition.h"

static int air_threshold(const u8*v,size_t n){ long h[256]={0},t=0; double s=0; long z=0;
  for(size_t i=0;i<n;i++){int x=v[i]; if(x){s+=x;z++;} if(x>=1&&x<=254){h[x]++;t++;}}
  if(t<256)return z?(int)(0.5*s/z+0.5):1; double sum=0; for(int i=1;i<=254;i++)sum+=(double)i*h[i];
  double sB=0,wB=0,b=-1; int thr=1; for(int k=1;k<=254;k++){wB+=h[k]; if(!wB)continue; double wF=t-wB; if(wF<=0)break;
    sB+=(double)k*h[k]; double mB=sB/wB,mF=(sum-sB)/wF,bb=wB*wF*(mB-mF)*(mB-mF); if(bb>b){b=bb;thr=k;}} return thr; }

/* Dominant eigenvector of a 3x3 symmetric matrix A=[xx,yy,zz,xy,xz,yz] via cyclic Jacobi.
 * Returns the eigenvector of the LARGEST eigenvalue (the across-sheet normal). */
static void eig3_dominant(const double A[6], double evec[3]){
  double a[3][3]={{A[0],A[3],A[4]},{A[3],A[1],A[5]},{A[4],A[5],A[2]}};
  double V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  for(int sweep=0;sweep<12;sweep++){
    double off=fabs(a[0][1])+fabs(a[0][2])+fabs(a[1][2]); if(off<1e-12)break;
    for(int p=0;p<2;p++)for(int q=p+1;q<3;q++){ if(fabs(a[p][q])<1e-15)continue;
      double th=0.5*(a[q][q]-a[p][p])/a[p][q];
      double t=(th>=0?1.0:-1.0)/(fabs(th)+sqrt(th*th+1.0));
      double c=1.0/sqrt(t*t+1.0),s=t*c;
      double app=a[p][p],aqq=a[q][q],apq=a[p][q];
      a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
      for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=a[k][p],akq=a[k][q];
        a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
      for(int k=0;k<3;k++){ double vkp=V[k][p],vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq; }
    }
  }
  int bi=0; double bl=a[0][0]; if(a[1][1]>bl){bl=a[1][1];bi=1;} if(a[2][2]>bl){bl=a[2][2];bi=2;}
  evec[0]=V[0][bi]; evec[1]=V[1][bi]; evec[2]=V[2][bi];
  double nl=sqrt(evec[0]*evec[0]+evec[1]*evec[1]+evec[2]*evec[2]); if(nl>1e-9){evec[0]/=nl;evec[1]/=nl;evec[2]/=nl;}
}

/* NaN-aware trilinear sample of a coarse winding volume (off-material = NaN; average finite corners) */
static double cw_trilin(const float*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0,acc=0,wt=0;
  int zs[2]={z0,z1},ys[2]={y0,y1},xs[2]={x0,x1}; double zw[2]={1-dz,dz},yw[2]={1-dy,dy},xw[2]={1-dx,dx};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){
    double vv=w[((size_t)zs[a]*ny+ys[b])*nx+xs[c]]; if(!isfinite(vv))continue;
    double gg=zw[a]*yw[b]*xw[c]; acc+=gg*vv; wt+=gg; }
  return wt>1e-9? acc/wt : NAN;
}

/* open-addressing pair hash: key (a<<32|b), a<b -> edge index, accumulate across/sheet/cnt */
typedef struct { u32 a,b; double sacross,ssheet; u32 cnt; } pedge;
int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [d_seed=4] [compact=20] [krepel=1] [kattract=1] [zc=mid] [cy=auto] [cx=auto] [pitch=8] [wgate=0.6]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int dseed=argc>10?atoi(argv[10]):4; float compact=argc>11?atof(argv[11]):20.0f;
  double krepel=argc>12?atof(argv[12]):1.0, kattract=argc>13?atof(argv[13]):1.0;
  int zc=argc>14?atoi(argv[14]):dz/2;
  // PHASE 1 winding gate (docs/touching-sheets-plan.md): merge supervoxels only if their winding
  // numbers agree. Winding is GEOMETRIC and non-circular: on an Archimedean spiral w=(r-a)/pitch, so
  // for two adjacent supervoxels dw = |r_i - r_j| / pitch (r = in-plane distance from umbilicus). A
  // tangential same-wrap neighbor has dr~0 (dw~0); a touching next wrap is ~pitch out (dw~1) even when
  // every LOCAL feature (across, sheetness) is identical -- the discriminator the affinity weight lacks.
  double cy=argc>15?atof(argv[15]):-1, cx=argc>16?atof(argv[16]):-1;   // umbilicus (region-local); auto if <0
  double pitch=argc>17?atof(argv[17]):8.0;                              // radial wrap spacing in voxels (L2~8)
  double wgate=argc>18?atof(argv[18]):0.6;                              // |dw|>=wgate -> different wrap (mutex)
  // Phase 1b: a COARSE WINDING FIELD (sheet_sep3d _vol.f32) gives dw that respects sheet DEFORMATION,
  // where raw geometric r/pitch fails (deformed wraps' radius-from-center swings around the wrap). This
  // is Thaumato/FASP's actual signal (winding = unwrapped angle the global solve provides). Pass
  // priorvol+priorlod to use it; else fall back to geometric r/pitch.
  const char*pvf=argc>19?argv[19]:""; int plod=argc>20?atoi(argv[20]):lod;
  float*cw=NULL; int cnz=0,cny=0,cnx=0; long cz0=0,cy0=0,cx0=0;
  if(pvf[0]){ FILE*pf=fopen(pvf,"rb");
    if(pf){ int hd[6]; if(fread(hd,sizeof(int),6,pf)==6){ cnz=hd[0];cny=hd[1];cnx=hd[2];cz0=hd[3];cy0=hd[4];cx0=hd[5];
      size_t cn=(size_t)cnz*cny*cnx; cw=malloc(cn*sizeof(float));
      if(fread(cw,sizeof(float),cn,pf)!=cn){ free(cw); cw=NULL; } } fclose(pf);
      if(cw) fprintf(stderr,"coarse winding %dx%dx%d @lod%d origin(%ld,%ld,%ld)\n",cnz,cny,cnx,plod,cz0,cy0,cx0); }
    if(!cw) fprintf(stderr,"WARN: could not read priorvol %s -- using geometric dw\n",pvf); }
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; double scy=0,scx=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr; nmat+=mask[p];
    if(mask[p]){ int x=p%dx,y=(p/dx)%dy; scy+=y; scx+=x; } }
  if(cy<0) cy = nmat? scy/nmat : dy/2.0;   // auto umbilicus = material centroid (≈ center for a centered crop)
  if(cx<0) cx = nmat? scx/nmat : dx/2.0;
  fprintf(stderr,"region %dx%dx%d air<%d, %ld material; umbilicus (y=%.1f x=%.1f) pitch=%.2f wgate=%.2f\n",
    dz,dy,dx,athr,nmat,cy,cx,pitch,wgate);

  // sheetness + per-voxel normals
  f32*sh=malloc(nn*sizeof(f32)),*nrm=malloc(3*nn*sizeof(f32)),*vf=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=1.5f;
  st_sheet_detect(vf,dz,dy,dx,&sp,sh,nrm);
  double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p]; if(smx<1e-6)smx=1;
  for(size_t p=0;p<nn;p++){ sh[p]=(f32)(sh[p]/smx); vf[p]=(f32)(sh[p]*255.0); }  // sh in [0,1]; vf=0..255 for SNIC

  // SNIC supervoxels on the sheetness field (they hug the sheets)
  int nsp=snic_superpixel_count(dz,dy,dx,dseed);
  u32*labels=calloc(nn,sizeof(u32)); Superpixel*sps=calloc((size_t)nsp+1,sizeof(Superpixel));
  int ov=snic(vf,dz,dy,dx,dseed,compact,80.0f,160.0f,labels,sps); free(vf);
  fprintf(stderr,"SNIC: d_seed=%d -> %d supervoxels (neigh-overflow=%d)\n",dseed,nsp,ov);

  // 1) per-supervoxel orientation: sheetness-weighted tensor sum Sum w*(n n^T), dominant eigvec = normal
  double*T=calloc((size_t)(nsp+1)*6,sizeof(double));
  for(size_t p=0;p<nn;p++){ if(!mask[p])continue; u32 L=labels[p]; if(!L)continue;
    double w=sh[p]; if(w<1e-3)continue; double*Tp=T+(size_t)L*6;
    double nx=nrm[3*p+0],ny=nrm[3*p+1],nz=nrm[3*p+2];
    Tp[0]+=w*nx*nx; Tp[1]+=w*ny*ny; Tp[2]+=w*nz*nz; Tp[3]+=w*nx*ny; Tp[4]+=w*nx*nz; Tp[5]+=w*ny*nz; }
  free(nrm);
  float*SN=calloc((size_t)(nsp+1)*3,sizeof(float)); long oriented=0;
  for(int L=1;L<=nsp;L++){ double*Tp=T+(size_t)L*6;
    if(Tp[0]+Tp[1]+Tp[2]<1e-9)continue; double ev[3]; eig3_dominant(Tp,ev);
    SN[3*L+0]=(float)ev[0]; SN[3*L+1]=(float)ev[1]; SN[3*L+2]=(float)ev[2]; oriented++; }
  free(T); fprintf(stderr,"oriented %ld/%d supervoxels (sheetness-weighted tensor)\n",oriented,nsp);

  // 2) supervoxel affinity: accumulate across/sheet over all boundary faces per adjacent pair
  size_t hcap=1; while(hcap < (size_t)nsp*16) hcap<<=1; if(hcap<1024)hcap=1024;
  s32*htab=malloc(hcap*sizeof(s32)); for(size_t i=0;i<hcap;i++) htab[i]=-1;
  int ecap=4096,ne=0; pedge*E=malloc((size_t)ecap*sizeof(pedge));
  #define HMASK (hcap-1)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
    if(!mask[p])continue; u32 Li=labels[p]; if(!Li)continue;
    for(int ax=0;ax<3;ax++){ int zz=z+(ax==2),yy=y+(ax==1),xx=x+(ax==0);
      if(xx>=dx||yy>=dy||zz>=dz)continue; size_t q=((size_t)zz*dy+yy)*dx+xx;
      if(!mask[q])continue; u32 Lj=labels[q]; if(!Lj||Lj==Li)continue;     // only inter-supervoxel faces
      u32 a=Li<Lj?Li:Lj, b=Li<Lj?Lj:Li;
      double across=0.5*(fabs((double)SN[3*Li+ax])+fabs((double)SN[3*Lj+ax]));  // contact along axis ax
      double sb=0.5*((double)sh[p]+(double)sh[q]);
      u64 key=((u64)a<<32)|b; size_t h=(size_t)((key*0x9E3779B97F4A7C15ull)>>40)&HMASK; int idx;
      for(;;){ idx=htab[h]; if(idx<0){ if(ne==ecap){ecap*=2; E=realloc(E,(size_t)ecap*sizeof(pedge));}
          idx=ne++; E[idx].a=a;E[idx].b=b;E[idx].sacross=0;E[idx].ssheet=0;E[idx].cnt=0; htab[h]=idx; break; }
        if(E[idx].a==a&&E[idx].b==b) break; h=(h+1)&HMASK; }
      E[idx].sacross+=across; E[idx].ssheet+=sb; E[idx].cnt++; }
  }
  free(htab);
  fprintf(stderr,"supervoxel RAG: %d adjacent pairs\n",ne);

  // per-supervoxel in-plane polar coords about the umbilicus (centroid sps[L].x/.y; x fastest, y=dy axis)
  float*R=calloc((size_t)(nsp+1),sizeof(float)),*TH=calloc((size_t)(nsp+1),sizeof(float));
  for(int L=1;L<=nsp;L++) if(sps[L].n){ double ry=sps[L].y-cy, rx=sps[L].x-cx;
    R[L]=(float)sqrt(ry*ry+rx*rx); TH[L]=(float)atan2(ry,rx); }
  // per-supervoxel winding sampled from the coarse field at the centroid (mapped lod->plod grid)
  float*W=NULL; long nwf=0; if(cw){ W=malloc((size_t)(nsp+1)*sizeof(float)); double scl=ldexp(1.0,lod-plod);
    for(int L=1;L<=nsp;L++){ W[L]=NAN; if(!sps[L].n)continue;
      double zc2=(z0+sps[L].z)*scl-cz0, yc2=(y0+sps[L].y)*scl-cy0, xc2=(x0+sps[L].x)*scl-cx0;
      double wv=cw_trilin(cw,cnz,cny,cnx,zc2,yc2,xc2); if(isfinite(wv)){ W[L]=(float)wv; nwf++; } }
    fprintf(stderr,"coarse winding sampled at %ld/%d supervoxel centroids\n",nwf,nsp); }
  // Archimedean winding w = r/pitch - s*theta/2pi (s=handedness): along ONE wrap r grows by a full pitch
  // per turn, so the theta term CANCELS the radial growth -> tangential same-wrap neighbors get dw~0. Pure
  // dr/pitch (no theta) was the bug: it cut every wrap radially as the spiral climbs. Pick s by minimizing
  // sum|dw| over edges (most adjacent pairs are tangential same-wrap, dw->0 only with the right handedness).
  #define DWE(e,s) (((double)R[E[e].a]-(double)R[E[e].b])/(pitch>1e-6?pitch:1) \
    - (s)*atan2(sin((double)TH[E[e].a]-TH[E[e].b]),cos((double)TH[E[e].a]-TH[E[e].b]))/(2*M_PI))
  double sp_=0,sm_=0; for(int e=0;e<ne;e++){ sp_+=fabs(DWE(e,1.0)); sm_+=fabs(DWE(e,-1.0)); }
  double hand = sp_<=sm_? 1.0 : -1.0;
  fprintf(stderr,"winding handedness s=%+.0f (sum|dw|: +%.0f / -%.0f)\n",hand,sp_,sm_);

  // build the signed graph with the PHASE 1 WINDING GATE. |dw|>=wgate => different wrap => hard MUTEX,
  // regardless of how high across*sheetness is (defeats a fused touch the local affinity can't see).
  // |dw|<wgate => same wrap => attract, weighted by (1-across) so MWS merges the tangentially-coherent
  // supervoxels of a wrap FIRST, then the winding mutexes separate the wraps.
  sgraph g; g.nnodes=nsp; g.nedges=ne; g.edges=malloc((size_t)ne*sizeof(sg_edge));
  long nrep=0, ngate=0, nold=0, nfield=0; double wsum=0;
  for(int e=0;e<ne;e++){ double across=E[e].sacross/E[e].cnt, sb=E[e].ssheet/E[e].cnt;
    double dw; u32 a=E[e].a,b=E[e].b;
    if(W && isfinite(W[a]) && isfinite(W[b])){ dw=fabs((double)W[a]-(double)W[b]); nfield++; }  // field (deformation-aware)
    else dw=fabs(DWE(e,hand));                                                                   // geometric fallback
    double w_old=kattract*(1.0-across) - krepel*(across*sb); if(w_old<0)nold++;   // pre-gate (diagnostic)
    double w;
    if(dw>=wgate){ w = -(krepel+0.5); ngate++; }       // winding says different wrap -> repel/mutex
    else         { w = kattract*(1.0-across); }         // same wrap -> attract (tangential-first)
    g.edges[e].a=a-1; g.edges[e].b=b-1; g.edges[e].w=(f32)w; if(w<0)nrep++; wsum+=w; }
  free(E); free(R); free(TH); if(cw)free(cw);
  fprintf(stderr,"dw source: %ld field / %ld geometric\n",nfield,ne-nfield);
  fprintf(stderr,"signed edges: %d (%ld repulsive: %ld by WINDING-GATE dw>=%.2f vs %ld by old across*sheet; mean w=%.3f)\n",
    ne,nrep,ngate,wgate,nold,ne?wsum/ne:0);

  // partition
  u32*seg=malloc((size_t)nsp*sizeof(u32));
  int ncl=mws_partition(&g,seg);
  // count non-trivial segments (>1 supervoxel) for a fragmentation read
  u32*scnt=calloc((size_t)ncl,sizeof(u32)); for(int L=1;L<=nsp;L++) if(sps[L].n) scnt[seg[L-1]]++;
  long big=0; for(int c=0;c<ncl;c++) if(scnt[c]>1) big++;
  fprintf(stderr,"Mutex-Watershed: %d segments over %d supervoxels (%ld multi-supervoxel)\n",ncl,nsp,big);
  free(scnt);

  // color zc slice by segment
  u8*rgb=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int gg=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)gg;
    if(mask[p] && labels[p]){ u32 S=seg[labels[p]-1],h=S*2654435761u;
      rgb[o]=(u8)(70+(h&150)); rgb[o+1]=(u8)(70+((h>>8)&150)); rgb[o+2]=(u8)(70+((h>>16)&150)); } }
  char fn[700]; snprintf(fn,sizeof fn,"%s_svseg.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }

  // DIAGNOSTIC: color zc slice by integer WINDING BAND floor(field winding) -- does the field ALONE give
  // clean concentric wraps (Family A: winding number IS the instance label), independent of MWS?
  if(W){ u8*rb=calloc((size_t)dy*dx*3,1);
    for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
      int gg=v[p]/4; rb[o]=rb[o+1]=rb[o+2]=(u8)gg;
      if(mask[p]&&labels[p]&&isfinite(W[labels[p]])){ int S=(int)floorf(W[labels[p]]);
        // ORDERED cycling palette (band%6): adjacent wraps differ but concentric continuity is visible
        // (unlike a random hash, which makes clean thin bands look like noise).
        static const u8 pal[6][3]={{230,60,60},{240,180,40},{60,200,80},{50,180,220},{90,90,230},{220,90,210}};
        int q=((S%6)+6)%6; rb[o]=pal[q][0]; rb[o+1]=pal[q][1]; rb[o+2]=pal[q][2]; } }
    char fb[700]; snprintf(fb,sizeof fb,"%s_band.ppm",outp); FILE*fbf=fopen(fb,"wb");
    if(fbf){ fprintf(fbf,"P6\n%d %d\n255\n",dx,dy); fwrite(rb,1,(size_t)dy*dx*3,fbf); fclose(fbf); fprintf(stderr,"wrote %s\n",fb); }
    free(rb); }
  free(rgb);free(seg);free(g.edges);free(SN);free(labels);free(sps);free(sh);free(mask);free(v); if(W)free(W); mca_close(arc); return 0;
}
