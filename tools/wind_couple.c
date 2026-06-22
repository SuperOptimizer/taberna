/* wind_couple.c -- z-COUPLED 3D winding solve over a stack of per-slice winding fields.
 *
 * sheet_sep solves the winding per z-plane independently; the planes are individually
 * good but don't agree across z (only ~59% of voxels within one wrap on adjacent L2
 * slices) because each plane anchors its winding origin independently and has its own
 * local switch defects. This tool couples them:
 *
 *   1. load N per-slice _wind.f32 dumps (hdr {d,z0,y0,x0}+{cyf,cxf,pitch}+d*d float,
 *      NAN off-material -- so barriers/air are simply absent and never couple).
 *   2. ALIGN the integer winding origin of each slice to its predecessor (the planes
 *      differ by whole wraps; round the median difference).
 *   3. 3D anisotropic red-black SOR: each material voxel relaxes toward its in-plane
 *      (x,y) AND z neighbours, plus a DATA TERM holding it near its validated 2D value.
 *      The z-coupling forces winding continuity through the volume; the data term keeps
 *      it from drifting off the 2D solution. Barriers/air (NAN) carry no edge -> the
 *      in-plane sheet separation is preserved exactly.
 *
 * Output: coherence stats (fraction within 1/2 and 1 wrap of the per-voxel z-median,
 * before align / after align / after couple), a coupled volume (OUT_z<z>.f32), and a
 * (z,radius) re-slice PPM (smooth vertical colour == z-coherent).
 *
 * usage: wind_couple OUTBASE f0.f32 f1.f32 ...   [--lam L] [--wz W] [--wxy W] [--it N]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef unsigned char u8;
typedef struct { int d,z0,y0,x0; double cyf,cxf,pitch; float *w; } slice;

static int cmpz(const void*a,const void*b){ const slice*x=a,*y=b; return x->z0-y->z0; }
static int cmpf(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return (x>y)-(x<y); }

static int load_slice(const char*path, slice*s){
  FILE*f=fopen(path,"rb"); if(!f){ fprintf(stderr,"open %s fail\n",path); return -1; }
  int hdr[4]; double geo[3];
  if(fread(hdr,sizeof(int),4,f)!=4 || fread(geo,sizeof(double),3,f)!=3){ fclose(f); return -1; }
  s->d=hdr[0]; s->z0=hdr[1]; s->y0=hdr[2]; s->x0=hdr[3];
  s->cyf=geo[0]; s->cxf=geo[1]; s->pitch=geo[2];
  size_t nn=(size_t)s->d*s->d; s->w=malloc(nn*sizeof(float));
  if(fread(s->w,sizeof(float),nn,f)!=nn){ fclose(f); free(s->w); return -1; }
  fclose(f); return 0;
}

// coherence: fraction of (slice,voxel) within tol of that voxel's z-median, over material
static void coherence(slice*S,int N,size_t nn,double*w05,double*w10){
  long in05=0,in10=0,tot=0;
  float *col=malloc(N*sizeof(float));
  for(size_t p=0;p<nn;p++){
    int m=0; for(int i=0;i<N;i++) if(!isnan(S[i].w[p])) col[m++]=S[i].w[p];
    if(m<2) continue;
    for(int a=0;a<m;a++)for(int b=a+1;b<m;b++) if(col[b]<col[a]){float t=col[a];col[a]=col[b];col[b]=t;}
    double med = (m&1)? col[m/2] : 0.5*(col[m/2-1]+col[m/2]);
    for(int i=0;i<N;i++){ float w=S[i].w[p]; if(isnan(w))continue; double dv=fabs(w-med); tot++; if(dv<0.5)in05++; if(dv<1.0)in10++; }
  }
  free(col);
  *w05=100.0*in05/(tot?tot:1); *w10=100.0*in10/(tot?tot:1);
}

int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s OUTBASE f0.f32 ... [--lam L][--wz W][--wxy W][--it N]\n",argv[0]); return 2; }
  const char*outb=argv[1];
  double LAM=0.5, WZ=1.0, WXY=0.25; int ITER=200;
  char*files[4096]; int N=0;
  for(int i=2;i<argc;i++){
    if(!strcmp(argv[i],"--lam")&&i+1<argc) LAM=atof(argv[++i]);
    else if(!strcmp(argv[i],"--wz")&&i+1<argc) WZ=atof(argv[++i]);
    else if(!strcmp(argv[i],"--wxy")&&i+1<argc) WXY=atof(argv[++i]);
    else if(!strcmp(argv[i],"--it")&&i+1<argc) ITER=atoi(argv[++i]);
    else if(N<4096) files[N++]=argv[i];
  }
  if(N<2){ fprintf(stderr,"need >=2 slices\n"); return 2; }
  slice *S=calloc(N,sizeof(slice)); int nok=0;
  for(int i=0;i<N;i++) if(load_slice(files[i],&S[nok])==0) nok++;
  N=nok; if(N<2){ fprintf(stderr,"loaded <2 slices\n"); return 1; }
  qsort(S,N,sizeof(slice),cmpz);
  int d=S[0].d; size_t nn=(size_t)d*d;
  for(int i=1;i<N;i++) if(S[i].d!=d){ fprintf(stderr,"slice d mismatch\n"); return 1; }
  fprintf(stderr,"%d slices d=%d  z=[%d..%d]  pitch=%.1f\n",N,d,S[0].z0,S[N-1].z0,S[0].pitch);

  double a05,a10; coherence(S,N,nn,&a05,&a10);
  printf("coherence BEFORE align: within0.5=%.0f%% within1=%.0f%%\n",a05,a10);

  // ALIGN integer winding origin: slice i -> match running predecessor (median diff).
  for(int i=1;i<N;i++){
    double sum=0; long c=0;
    // collect diffs; median via cheap approach: accumulate into a small histogram of rounded diffs
    long hist[401]={0}; // diff in [-200,200]
    for(size_t p=0;p<nn;p++){ float a=S[i-1].w[p],b=S[i].w[p]; if(isnan(a)||isnan(b))continue;
      long dd=lround(a-b); if(dd<-200)dd=-200; if(dd>200)dd=200; hist[dd+200]++; c++; (void)sum; }
    if(c<1000) continue;
    long best=0,bc=-1; for(int h=0;h<401;h++) if(hist[h]>bc){bc=hist[h];best=h-200;}  // mode of integer diff
    for(size_t p=0;p<nn;p++) if(!isnan(S[i].w[p])) S[i].w[p]+=best;
  }
  coherence(S,N,nn,&a05,&a10);
  printf("coherence AFTER  align: within0.5=%.0f%% within1=%.0f%%\n",a05,a10);

  // ROBUST z-median pre-vote (window +/-2): the per-slice defects are INTEGER-JUMP
  // outliers (a voxel mislabelled by whole wraps). Linear diffusion AVERAGES them
  // (wrong); the median REJECTS them. So vote first -> robust target `orig`, then couple
  // toward it (the coupling adds spatial propagation the pure median lacks).
  float **orig=malloc(N*sizeof(float*));
  for(int i=0;i<N;i++) orig[i]=malloc(nn*sizeof(float));
  #pragma omp parallel
  { float col[64];
    #pragma omp for schedule(static)
    for(int i=0;i<N;i++){ int lo=i-2<0?0:i-2, hi=i+2>=N?N-1:i+2;
      for(size_t p=0;p<nn;p++){ if(isnan(S[i].w[p])){ orig[i][p]=NAN; continue; }
        int m=0; for(int j=lo;j<=hi;j++){ float w=S[j].w[p]; if(!isnan(w))col[m++]=w; }
        for(int a=0;a<m;a++)for(int b=a+1;b<m;b++) if(col[b]<col[a]){float t=col[a];col[a]=col[b];col[b]=t;}
        orig[i][p]= (m&1)? col[m/2] : 0.5f*(col[m/2-1]+col[m/2]); } } }
  for(int i=0;i<N;i++) memcpy(S[i].w,orig[i],nn*sizeof(float));   // robust initialisation
  coherence(S,N,nn,&a05,&a10);
  printf("coherence AFTER  vote : within0.5=%.0f%% within1=%.0f%%\n",a05,a10);

  // 3D anisotropic red-black SOR. Material = !isnan(orig). Barriers/air are NAN -> no edge.
  const double omega=1.6;
  for(int it=0;it<ITER;it++) for(int color=0;color<2;color++){
    #pragma omp parallel for schedule(static)
    for(int i=0;i<N;i++){
      float *w=S[i].w, *o=orig[i];
      float *wm=(i>0)?S[i-1].w:NULL, *wp=(i<N-1)?S[i+1].w:NULL;
      for(int y=0;y<d;y++)for(int x=0;x<d;x++){ if(((x+y+i)&1)!=color)continue;
        size_t p=(size_t)y*d+x; if(isnan(o[p]))continue;
        double wp0=w[p], s=0,ws=0,nv;
        // ROBUST (Tukey) edge weights: a neighbour far from the current estimate (an
        // integer-jump outlier) is down-weighted, so it cannot drag this voxel off-wrap.
        #define ROB(nb,wt) do{ double rr=((nb)-wp0)/0.6; double rw=(wt)/(1.0+rr*rr); s+=rw*(nb); ws+=rw; }while(0)
        if(x>0   && !isnan(w[p-1])){ nv=w[p-1]; ROB(nv,WXY); }
        if(x<d-1 && !isnan(w[p+1])){ nv=w[p+1]; ROB(nv,WXY); }
        if(y>0   && !isnan(w[p-d])){ nv=w[p-d]; ROB(nv,WXY); }
        if(y<d-1 && !isnan(w[p+d])){ nv=w[p+d]; ROB(nv,WXY); }
        if(wm && !isnan(wm[p])){ nv=wm[p]; ROB(nv,WZ); }
        if(wp && !isnan(wp[p])){ nv=wp[p]; ROB(nv,WZ); }
        #undef ROB
        s+=LAM*o[p]; ws+=LAM;                               // data term: stay near robust (voted) target
        if(ws>1e-9){ double tgt=s/ws; w[p]=(float)(w[p]+omega*(tgt-w[p])); }
      }
    }
  }
  coherence(S,N,nn,&a05,&a10);
  printf("coherence AFTER  couple: within0.5=%.0f%% within1=%.0f%%  (lam=%.2f wz=%.2f wxy=%.2f it=%d)\n",a05,a10,LAM,WZ,WXY,ITER);

  // write coupled volume + a (z,radius) re-slice PPM at a fixed angle
  for(int i=0;i<N;i++){ char fn[700]; snprintf(fn,sizeof fn,"%s_z%d.f32",outb,S[i].z0);
    FILE*f=fopen(fn,"wb"); if(f){ int hdr[4]={d,S[i].z0,S[i].y0,S[i].x0}; double geo[3]={S[i].cyf,S[i].cxf,S[i].pitch};
      fwrite(hdr,sizeof(int),4,f); fwrite(geo,sizeof(double),3,f); fwrite(S[i].w,sizeof(float),nn,f); fclose(f); } }
  { double ang=0.6, uy=sin(ang), ux=cos(ang); double cyf=S[0].cyf,cxf=S[0].cxf; int R=(int)(0.45*d);
    // ROBUST ramp range: sample winding values, take 2/98 percentiles (raw min/max are
    // blown out by a few integer-jump outliers -> everything else maps to one colour).
    int cap=200000; float*samp=malloc(cap*sizeof(float)); int ns=0;
    for(int i=0;i<N && ns<cap;i++)for(size_t p=0;p<nn && ns<cap;p+=37){ float w=S[i].w[p]; if(!isnan(w))samp[ns++]=w; }
    qsort(samp,ns,sizeof(float),cmpf);
    double lo= ns?samp[(int)(0.02*ns)]:0, hi= ns?samp[(int)(0.98*ns)]:1; free(samp);
    double sp=(hi>lo)?hi-lo:1.0;
    u8*img=calloc((size_t)N*R*3,1);
    for(int i=0;i<N;i++)for(int k=0;k<R;k++){ int yy=(int)lround(cyf+uy*k),xx=(int)lround(cxf+ux*k);
      if(yy<0||yy>=d||xx<0||xx>=d)continue; float w=S[i].w[(size_t)yy*d+xx]; if(isnan(w))continue;
      double t=(w-lo)/sp; if(t<0)t=0;if(t>1)t=1; size_t q=((size_t)i*R+k)*3; img[q]=(u8)(255*t); img[q+2]=(u8)(255*(1-t)); }
    char fn[700]; snprintf(fn,sizeof fn,"%s_reslice.ppm",outb); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P6\n%d %d\n255\n",R,N); fwrite(img,1,(size_t)N*R*3,f); fclose(f); }
    free(img); printf("wrote %s_reslice.ppm (rows=z cols=radius; smooth vertical colour = coherent) and %d coupled f32\n",outb,N); }
  return 0;
}
