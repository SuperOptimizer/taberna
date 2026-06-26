/* wind_tifxyz — export the winding field as a VC3D tifxyz QuadSurface ("plug ours into VC3D"). The
 * coherent unrolled papyrus surface, parameterized by (u = spiral render-coord, v = z), each pixel
 * holding the source VOLUME coordinate (z,y,x) it came from. This is the global-winding answer to the
 * surface grower: no local front to fragment — every wrap is laid out by the smooth field.
 *
 * Same spiral render coord as unroll_wind: rc = w + atan2(y-cy,x-cx)/2pi. Scatter each material voxel
 * into its (z, u=(rc-wmin)*SAMP) pixel, MEAN-accumulate its (z0+z,y0+y,x0+x); -1 = unfilled. Writes
 * OUTDIR/{x,y,z}.tif (f32) + mask.tif (u8) + meta.json. Coords are in the archive's `lod` grid.
 *
 *   wind_tifxyz ARC lod WINDING.f32 OUTDIR [samp=160]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "annotate/umbilicus.h"
#include "tiff.h"

#define TWO_PI 6.283185307179586

static int auto_center(const char*arc,int lod,int z0,int y0,int x0,int dz,double*cy,double*cx){
  mca_handle*h=mca_open(arc); if(!h)return -1; int fz,fy,fx,nl; float ql; mca_handle_dims(h,&fz,&fy,&fx,&ql,&nl);
  int cl=5; double cs=(double)(1<<cl); int ccz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*c=mca_read(h,cl,0,0,0,ccz,ccy,ccx); if(!c){mca_close(h);return -1;}
  size_t n=(size_t)ccz*ccy*ccx; u8*m=malloc(n); for(size_t i=0;i<n;i++)m[i]=c[i]!=0; free(c);
  umbilicus umb; int rc=umbilicus_estimate(m,ccz,ccy,ccx,9,&umb); free(m); mca_close(h);
  if(rc)return -1; double ls=(double)(1<<lod); double coarse_z=(z0+dz/2.0)*ls/cs;
  f32 ucy,ucx; umbilicus_center(&umb,(f32)coarse_z,&ucy,&ucx); umbilicus_free(&umb);
  *cy=ucy*cs/ls-y0; *cx=ucx*cs/ls-x0; return 0;
}
static int cmp_f(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return x<y?-1:x>y?1:0; }
static int wtif(const char*dir,const char*name,const float*d,int W,int H){
  char p[1200]; snprintf(p,sizeof p,"%s/%s.tif",dir,name);
  tiff_volume v={.width=(unsigned)W,.height=(unsigned)H,.depth=1,.channels=1,.type=TIFF_F32,.data=(void*)d};
  return tiff_write_volume(p,&v);
}

int main(int argc,char**argv){
  if(argc<5){ fprintf(stderr,"usage: %s ARC lod WINDING.f32 OUTDIR [samp=160]\n",argv[0]); return 2; }
  const char*arc=argv[1]; int lod=atoi(argv[2]); const char*wf=argv[3],*outdir=argv[4];
  double SAMP=argc>5?atof(argv[5]):160.0;
  FILE*f=fopen(wf,"rb"); if(!f){fprintf(stderr,"open winding fail\n");return 1;}
  int hd[6]; if(fread(hd,sizeof(int),6,f)!=6){return 1;} int dz=hd[0],dy=hd[1],dx=hd[2],z0=hd[3],y0=hd[4],x0=hd[5];
  size_t nn=(size_t)dz*dy*dx; float*w=malloc(nn*sizeof(float)); if(fread(w,sizeof(float),nn,f)!=nn){return 1;} fclose(f);

  double cy,cx; if(auto_center(arc,lod,z0,y0,x0,dz,&cy,&cx)){fprintf(stderr,"auto-center fail\n");return 1;}
  fprintf(stderr,"umbilicus (local) = (%.1f,%.1f); winding %dx%dx%d @ lod%d (z%d y%d x%d)\n",cy,cx,dz,dy,dx,lod,z0,y0,x0);

  // render-coord range over finite voxels (1..99.5 pctile, like unroll_wind)
  double wmin,wmax; { float*t=malloc(nn*sizeof(float)); size_t k=0;
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
      if(!isfinite(w[p]))continue; double rc=w[p]+atan2((double)y-cy,(double)x-cx)/TWO_PI; t[k++]=(float)rc; }
    if(!k){fprintf(stderr,"no winding\n");return 1;} qsort(t,k,sizeof(float),cmp_f);
    wmin=t[(size_t)(0.010*k)]; wmax=t[(size_t)(0.995*(k-1))]; free(t); }
  int UW=(int)((wmax-wmin)*SAMP)+1; if(UW>120000)UW=120000;
  fprintf(stderr,"render-coord %.1f..%.1f (%.0f wraps) -> %d cols x %d rows\n",wmin,wmax,wmax-wmin,UW,dz);

  // each (z,u) bin collects a small sheet patch (~radial thickness x azimuthal span); a MEAN coord drifts
  // off the wavy sheet -> noise. Keep the RIDGE voxel (max CT = the sheet centerline) as the surface point.
  mca_handle*h=mca_open(arc); if(!h){fprintf(stderr,"open arc fail\n");return 1;}
  u8*ct=mca_read(h,lod,z0,y0,x0,dz,dy,dx); mca_close(h); if(!ct){fprintf(stderr,"CT read fail\n");return 1;}
  // EXTRACT mode (WIND_EXTRACT env): "ridge" keeps the single max-CT voxel per (z,u) bin -- sharp but
  // geometrically NOISY (the brightest voxel hops radially between columns -> a jumpy 3D manifold).
  // "centroid" (default) takes the CT-weighted mean position over the bin: the bin is thin in render-
  // coord (~one sheet), so the centroid stays on the sheet but is geometrically SMOOTH -- a clean
  // manifold like a VC3D grown surface, without the front-fragmentation holes. WIND_SMOOTH=N adds an
  // NxN median filter over the extracted (z,u) coords to remove residual single-pixel jumps.
  int centroid=1; { const char*e=getenv("WIND_EXTRACT"); if(e&&!strcmp(e,"ridge"))centroid=0; }
  int smooth=0; { const char*e=getenv("WIND_SMOOTH"); if(e)smooth=atoi(e); if(smooth<0)smooth=0; }
  size_t no=(size_t)dz*UW; float*X=malloc(no*sizeof(float)),*Y=malloc(no*sizeof(float)),*Z=malloc(no*sizeof(float)),*M=malloc(no*sizeof(float));
  double*ax=NULL,*ay=NULL,*az=NULL,*aw=NULL; float*best=NULL;
  if(centroid){ ax=calloc(no,sizeof(double));ay=calloc(no,sizeof(double));az=calloc(no,sizeof(double));aw=calloc(no,sizeof(double)); }
  else best=malloc(no*sizeof(float));
  for(size_t i=0;i<no;i++){ X[i]=Y[i]=Z[i]=-1; M[i]=0; if(best)best[i]=-1; }
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
    float wv=w[p]; if(!isfinite(wv))continue;
    double rc=wv+atan2((double)y-cy,(double)x-cx)/TWO_PI; int u=(int)((rc-wmin)*SAMP); if(u<0||u>=UW)continue;
    size_t o=(size_t)z*UW+u;
    if(centroid){ double cw=ct[p]; if(cw<=0)continue; ax[o]+=cw*(x0+x); ay[o]+=cw*(y0+y); az[o]+=cw*(z0+z); aw[o]+=cw; M[o]=1; }
    else if((float)ct[p]>best[o]){ best[o]=ct[p]; X[o]=x0+x; Y[o]=y0+y; Z[o]=z0+z; M[o]=1; } }
  if(centroid){ for(size_t o=0;o<no;o++) if(aw[o]>0){ X[o]=(float)(ax[o]/aw[o]); Y[o]=(float)(ay[o]/aw[o]); Z[o]=(float)(az[o]/aw[o]); }
    free(ax);free(ay);free(az);free(aw); }
  if(smooth>0){ // separable-ish NxN median over valid coords, per channel
    int r=smooth/2; float*tmp=malloc(no*sizeof(float)); float*win=malloc((size_t)(smooth*smooth)*sizeof(float));
    float*chs[3]={X,Y,Z};
    for(int c=0;c<3;c++){ float*S=chs[c];
      for(int z=0;z<dz;z++)for(int u=0;u<UW;u++){ size_t o=(size_t)z*UW+u; if(M[o]<=0){tmp[o]=S[o];continue;}
        int k=0; for(int dz2=-r;dz2<=r;dz2++)for(int du=-r;du<=r;du++){ int zz=z+dz2,uu=u+du; if(zz<0||uu<0||zz>=dz||uu>=UW)continue;
          size_t oo=(size_t)zz*UW+uu; if(M[oo]>0)win[k++]=S[oo]; }
        if(!k){tmp[o]=S[o];continue;} for(int a=1;a<k;a++){float v=win[a];int b=a-1;while(b>=0&&win[b]>v){win[b+1]=win[b];b--;}win[b+1]=v;}
        tmp[o]=win[k/2]; }
      for(size_t o=0;o<no;o++) S[o]=tmp[o]; }
    free(tmp);free(win); }
  // WIND_SNAP=R: hybrid sharpen. The smooth centroid gives a clean manifold POSITION but blurs CT
  // across the sheet thickness. Snap each point along the radial (umbilicus-outward) normal within
  // +/-R lod-voxels to the brightest CT voxel -- restores the sharp sheet centerline without
  // reintroducing the max-ridge's lateral jumpiness (the search is confined to the normal line).
  int snap=0; { const char*e=getenv("WIND_SNAP"); if(e)snap=atoi(e); if(snap<0)snap=0; }
  if(snap>0){ long moved=0;
    for(size_t o=0;o<no;o++){ if(M[o]<=0)continue;
      double yi=Y[o]-y0, xi=X[o]-x0; int zi=(int)lround(Z[o]-z0);
      double ry=yi-cy, rx=xi-cx, rn=sqrt(ry*ry+rx*rx); if(rn<1e-6||zi<0||zi>=dz)continue; ry/=rn; rx/=rn;
      double bestv=-1, bt=0;
      for(double t=-snap;t<=snap;t+=0.5){ int yy=(int)lround(yi+t*ry), xx=(int)lround(xi+t*rx);
        if(yy<0||xx<0||yy>=dy||xx>=dx)continue; double cv=ct[((size_t)zi*dy+yy)*dx+xx]; if(cv>bestv){bestv=cv;bt=t;} }
      if(bt!=0){ Y[o]=(float)(y0+yi+bt*ry); X[o]=(float)(x0+xi+bt*rx); moved++; } }
    fprintf(stderr,"snap R=%d: moved %ld points to radial CT ridge\n",snap,moved); }
  free(ct);
  long filled=0; for(size_t i=0;i<no;i++) filled+=(M[i]>0);
  fprintf(stderr,"extract=%s smooth=%d  filled %ld/%zu pixels (%.1f%%)\n",centroid?"centroid":"ridge",smooth,filled,no,100.0*filled/no);

  if(wtif(outdir,"x",X,UW,dz)||wtif(outdir,"y",Y,UW,dz)||wtif(outdir,"z",Z,UW,dz)||wtif(outdir,"mask",M,UW,dz)){
    fprintf(stderr,"tif write fail (does %s exist?)\n",outdir); return 1; }
  char mp[1200]; snprintf(mp,sizeof mp,"%s/meta.json",outdir); FILE*mf=fopen(mp,"w");
  if(mf){ fprintf(mf,"{\n  \"type\":\"seg\", \"source\":\"taberna wind_tifxyz\",\n  \"width\":%d, \"height\":%d, \"lod\":%d,\n  \"region\":[%d,%d,%d, %d,%d,%d],\n  \"render_coord\":[%.3f,%.3f]\n}\n",
      UW,dz,lod,z0,y0,x0,dz,dy,dx,wmin,wmax); fclose(mf); }
  fprintf(stderr,"wrote %s/{x,y,z,mask}.tif + meta.json (%dx%d)\n",outdir,UW,dz);
  free(w);if(best)free(best);free(X);free(Y);free(Z);free(M); return 0;
}
