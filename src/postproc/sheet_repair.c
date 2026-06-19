/* sheet_repair.c — see sheet_repair.h. */
#include "postproc/sheet_repair.h"
#include "eval/topo.h"   // cc_label

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))

/* Jacobi eigendecomposition of a 3x3 symmetric matrix. eval ascending, evec cols. */
static void eig3(double A[3][3], double eval[3], double evec[3][3]) {
  double a[3][3]; memcpy(a, A, sizeof a);
  double v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  for (int iter=0; iter<50; iter++) {
    // largest off-diagonal
    int p=0,q=1; double mx=fabs(a[0][1]);
    if (fabs(a[0][2])>mx){mx=fabs(a[0][2]);p=0;q=2;}
    if (fabs(a[1][2])>mx){mx=fabs(a[1][2]);p=1;q=2;}
    if (mx < 1e-12) break;
    double app=a[p][p],aqq=a[q][q],apq=a[p][q];
    double phi=0.5*atan2(2*apq, aqq-app);
    double c=cos(phi), s=sin(phi);
    for (int k=0;k<3;k++){ double akp=a[k][p],akq=a[k][q]; a[k][p]=c*akp - s*akq; a[k][q]=s*akp + c*akq; }
    for (int k=0;k<3;k++){ double apk=a[p][k],aqk=a[q][k]; a[p][k]=c*apk - s*aqk; a[q][k]=s*apk + c*aqk; }
    for (int k=0;k<3;k++){ double vkp=v[k][p],vkq=v[k][q]; v[k][p]=c*vkp - s*vkq; v[k][q]=s*vkp + c*vkq; }
  }
  for (int i=0;i<3;i++) eval[i]=a[i][i];
  int ord[3]={0,1,2};
  for (int i=0;i<3;i++) for (int j=i+1;j<3;j++) if (eval[ord[j]]<eval[ord[i]]){int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  double ev[3]; double evc[3][3];
  for (int i=0;i<3;i++){ ev[i]=eval[ord[i]]; for(int k=0;k<3;k++) evc[k][i]=v[k][ord[i]]; }
  memcpy(eval,ev,sizeof ev); memcpy(evec,evc,sizeof evc);
}

static void set_vox(u8 *out, int nz,int ny,int nx, int z,int y,int x) {
  size_t nynx=(size_t)ny*nx;
  if (z>=0&&z<nz&&y>=0&&y<ny&&x>=0&&x<nx) out[IDX(z,y,x)]=1;
}
/* 3D line raster (DDA) between two voxel centers — keeps the surface connected. */
static void draw_line(u8 *out, int nz,int ny,int nx, double z0,double y0,double x0, double z1,double y1,double x1) {
  double dz=z1-z0, dy=y1-y0, dx=x1-x0;
  double L=fabs(dz); if(fabs(dy)>L)L=fabs(dy); if(fabs(dx)>L)L=fabs(dx);
  int n=(int)L+1;
  for (int i=0;i<=n;i++){ double t=n?(double)i/n:0;
    set_vox(out,nz,ny,nx,(int)lround(z0+t*dz),(int)lround(y0+t*dy),(int)lround(x0+t*dx)); }
}

void sheet_repair(const u8 *in, u8 *out, int nz, int ny, int nx, int min_voxels) {
  size_t nynx=(size_t)ny*nx, N=(size_t)nz*nynx;
  memset(out, 0, N);
  u32 *lab=malloc(N*sizeof(u32));
  u32 nc=cc_label(in, nz,ny,nx, TOPO_CONN26, lab);
  if (!nc){ free(lab); return; }
  size_t *sz=calloc(nc+1,sizeof(size_t));
  for (size_t i=0;i<N;i++) sz[lab[i]]++;

  // gather voxel indices per component: pts[off[c] .. off[c+1]) for component c
  size_t *off=calloc(nc+2,sizeof(size_t));
  for (u32 c=1;c<=nc;c++) off[c+1]=off[c]+ (sz[c]>=(size_t)min_voxels? sz[c]:0);
  size_t *pts=malloc((off[nc+1]+1)*sizeof(size_t));
  size_t *cur=malloc((nc+2)*sizeof(size_t)); memcpy(cur,off,(nc+2)*sizeof(size_t));
  for (size_t i=0;i<N;i++){ u32 c=lab[i]; if(c && sz[c]>=(size_t)min_voxels) pts[cur[c]++]=i; }
  free(cur);

  for (u32 c=1;c<=nc;c++) {
    size_t a=off[c], b=off[c+1];
    size_t m=b-a;
    if (m < (size_t)min_voxels) continue;

    // --- PCA ---
    double mean[3]={0,0,0};
    for (size_t i=a;i<b;i++){ size_t v=pts[i]; int z=v/nynx,r=v%nynx,y=r/nx,x=r%nx; mean[0]+=z;mean[1]+=y;mean[2]+=x; }
    mean[0]/=m; mean[1]/=m; mean[2]/=m;
    double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (size_t i=a;i<b;i++){ size_t v=pts[i]; int z=v/nynx,r=v%nynx,y=r/nx,x=r%nx;
      double d[3]={z-mean[0],y-mean[1],x-mean[2]};
      for(int p=0;p<3;p++)for(int q=0;q<3;q++) cov[p][q]+=d[p]*d[q]; }
    double ev[3],evec[3][3]; eig3(cov,ev,evec);
    // normal = smallest eigenvalue dir (col 0); tangents = cols 1,2
    double nrm[3]={evec[0][0],evec[1][0],evec[2][0]};
    double t1[3]={evec[0][1],evec[1][1],evec[2][1]};
    double t2[3]={evec[0][2],evec[1][2],evec[2][2]};

    // --- project to (u,v,w) ---
    double umin=1e18,umax=-1e18,vmin=1e18,vmax=-1e18;
    double *U=malloc(m*sizeof(double)),*Vv=malloc(m*sizeof(double)),*W=malloc(m*sizeof(double));
    for (size_t i=a;i<b;i++){ size_t v=pts[i]; int z=v/nynx,r=v%nynx,y=r/nx,x=r%nx;
      double d[3]={z-mean[0],y-mean[1],x-mean[2]};
      double uu=d[0]*t1[0]+d[1]*t1[1]+d[2]*t1[2];
      double vv=d[0]*t2[0]+d[1]*t2[1]+d[2]*t2[2];
      double ww=d[0]*nrm[0]+d[1]*nrm[1]+d[2]*nrm[2];
      size_t j=i-a; U[j]=uu; Vv[j]=vv; W[j]=ww;
      if(uu<umin)umin=uu; if(uu>umax)umax=uu; if(vv<vmin)vmin=vv; if(vv>vmax)vmax=vv; }

    int pad=2;
    int nu=(int)(umax-umin)+1+2*pad, nv=(int)(vmax-vmin)+1+2*pad;
    if (nu<3||nv<3||(size_t)nu*nv>200000000ull){ free(U);free(Vv);free(W); continue; }
    size_t G=(size_t)nu*nv;
    double *h=calloc(G,sizeof(double)); int *cnt=calloc(G,sizeof(int));
    for (size_t j=0;j<m;j++){ int gi=(int)lround(U[j]-umin)+pad, gj=(int)lround(Vv[j]-vmin)+pad;
      size_t g=(size_t)gi*nv+gj; h[g]+=W[j]; cnt[g]++; }
    u8 *known=calloc(G,1);
    for (size_t g=0;g<G;g++) if(cnt[g]){ h[g]/=cnt[g]; known[g]=1; }

    // --- domain = morphological closing of `known` (box radius R): fills porosity
    //     gaps up to ~2R as interior while keeping the sheet outline. This is the
    //     key fix vs flood-from-border, which wrongly treated border-connected
    //     porosity as "outside". ---
    int R=6;
    u8 *dil=calloc(G,1), *dom=calloc(G,1);
    for (int gi=0;gi<nu;gi++) for(int gj=0;gj<nv;gj++){ if(!known[(size_t)gi*nv+gj]) continue;
      for(int di=-R;di<=R;di++){int ii=gi+di; if(ii<0||ii>=nu)continue;
        for(int dj=-R;dj<=R;dj++){int jj=gj+dj; if(jj<0||jj>=nv)continue; dil[(size_t)ii*nv+jj]=1; }}}
    for (int gi=0;gi<nu;gi++) for(int gj=0;gj<nv;gj++){ // erode dil -> dom
      int keep=1; for(int di=-R;di<=R&&keep;di++){int ii=gi+di; if(ii<0||ii>=nu){keep=0;break;}
        for(int dj=-R;dj<=R;dj++){int jj=gj+dj; if(jj<0||jj>=nv||!dil[(size_t)ii*nv+jj]){keep=0;break;}}}
      dom[(size_t)gi*nv+gj]=keep; }
    free(dil);

    // --- Laplace inpaint interior gaps (Gauss-Seidel) ---
    for (int it=0; it<300; it++){ double md=0;
      for (int gi=1;gi<nu-1;gi++) for(int gj=1;gj<nv-1;gj++){ size_t g=(size_t)gi*nv+gj;
        if(known[g]||!dom[g]) continue;
        double nv4=(h[g-nv]+h[g+nv]+h[g-1]+h[g+1])*0.25; double dd=fabs(nv4-h[g]); if(dd>md)md=dd; h[g]=nv4; }
      if (md<1e-3) break;
    }

    // --- rasterize watertight: draw each domain cell, connect to +u,+v neighbors ---
    for (int gi=0;gi<nu;gi++) for(int gj=0;gj<nv;gj++){ size_t g=(size_t)gi*nv+gj; if(!dom[g])continue;
      double uu=(gi-pad)+umin, vv=(gj-pad)+vmin, ww=h[g];
      double wz=mean[0]+uu*t1[0]+vv*t2[0]+ww*nrm[0];
      double wy=mean[1]+uu*t1[1]+vv*t2[1]+ww*nrm[1];
      double wx=mean[2]+uu*t1[2]+vv*t2[2]+ww*nrm[2];
      set_vox(out,nz,ny,nx,(int)lround(wz),(int)lround(wy),(int)lround(wx));
      // connect to right/bottom domain neighbors to seal the surface
      if (gi+1<nu){ size_t g2=g+nv; if(dom[g2]){ double u2=(gi+1-pad)+umin,w2=h[g2];
        draw_line(out,nz,ny,nx, wz,wy,wx,
          mean[0]+u2*t1[0]+vv*t2[0]+w2*nrm[0], mean[1]+u2*t1[1]+vv*t2[1]+w2*nrm[1], mean[2]+u2*t1[2]+vv*t2[2]+w2*nrm[2]); } }
      if (gj+1<nv){ size_t g2=g+1; if(dom[g2]){ double v2=(gj+1-pad)+vmin,w2=h[g2];
        draw_line(out,nz,ny,nx, wz,wy,wx,
          mean[0]+uu*t1[0]+v2*t2[0]+w2*nrm[0], mean[1]+uu*t1[1]+v2*t2[1]+w2*nrm[1], mean[2]+uu*t1[2]+v2*t2[2]+w2*nrm[2]); } }
    }
    free(h);free(cnt);free(known);free(dom);free(U);free(Vv);free(W);
  }
  free(lab);free(sz);free(off);free(pts);
}

void sheet_repair_windowed(const u8 *in, u8 *out, int nz, int ny, int nx,
                           int win, int overlap, int min_voxels) {
  size_t nynx=(size_t)ny*nx, N=(size_t)nz*nynx;
  memset(out, 0, N);
  int stride = win - overlap; if (stride < 1) stride = 1;
  u8 *si=malloc((size_t)win*win*win), *so=malloc((size_t)win*win*win);
  for (int z0=0; z0<nz; z0+=stride)
  for (int y0=0; y0<ny; y0+=stride)
  for (int x0=0; x0<nx; x0+=stride) {
    int dz=win<nz-z0?win:nz-z0, dy=win<ny-y0?win:ny-y0, dx=win<nx-x0?win:nx-x0;
    if (dz<8||dy<8||dx<8) continue;
    size_t any=0;
    for (int z=0;z<dz;z++) for(int y=0;y<dy;y++) for(int x=0;x<dx;x++){
      u8 v=in[(size_t)(z0+z)*nynx+(size_t)(y0+y)*nx+(x0+x)];
      si[((size_t)z*dy+y)*dx+x]=v; any+=v; }
    if (any < (size_t)min_voxels) continue;
    sheet_repair(si, so, dz,dy,dx, min_voxels);
    for (int z=0;z<dz;z++) for(int y=0;y<dy;y++) for(int x=0;x<dx;x++)
      if (so[((size_t)z*dy+y)*dx+x]) out[(size_t)(z0+z)*nynx+(size_t)(y0+y)*nx+(x0+x)]=1;
  }
  free(si);free(so);
}
