/* hessian.c — see hessian.h. Multi-scale Frangi plate detector. */
#include "segmentation/hessian.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))
#define CL(v,n) ((v)<0?0:((v)>=(n)?(n)-1:(v)))

/* separable Gaussian blur (reflect borders) */
static void gblur(const f32 *src, f32 *dst, int nz, int ny, int nx, double sigma) {
  size_t nynx=(size_t)ny*nx, n=(size_t)nz*nynx;
  int r=(int)(3.0*sigma+0.5); if(r<1)r=1;
  double *k=malloc((2*r+1)*sizeof(double)); double s=0;
  for(int t=-r;t<=r;t++){ k[t+r]=exp(-0.5*(double)t*t/(sigma*sigma)); s+=k[t+r]; }
  for(int t=0;t<2*r+1;t++) k[t]/=s;
  f32 *a=malloc(n*sizeof(f32)), *b=malloc(n*sizeof(f32));
  memcpy(a,src,n*sizeof(f32));
  // x
  #pragma omp parallel for schedule(static)
  for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){ double acc=0;
    for(int t=-r;t<=r;t++) acc+=k[t+r]*a[IDX(z,y,CL(x+t,nx))]; b[IDX(z,y,x)]=(f32)acc; }
  // y
  #pragma omp parallel for schedule(static)
  for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){ double acc=0;
    for(int t=-r;t<=r;t++) acc+=k[t+r]*b[IDX(z,CL(y+t,ny),x)]; a[IDX(z,y,x)]=(f32)acc; }
  // z
  #pragma omp parallel for schedule(static)
  for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){ double acc=0;
    for(int t=-r;t<=r;t++) acc+=k[t+r]*a[IDX(CL(z+t,nz),y,x)]; dst[IDX(z,y,x)]=(f32)acc; }
  free(a);free(b);free(k);
}

/* eigen-decomposition of symmetric 3x3 (Jacobi). eval/evec by |eval| ascending. */
static void eig3m(double A[3][3], double ev[3], double V[3][3]){
  double a[3][3]; memcpy(a,A,sizeof a);
  double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  for(int it=0;it<50;it++){
    int p=0,q=1; double mx=fabs(a[0][1]);
    if(fabs(a[0][2])>mx){mx=fabs(a[0][2]);p=0;q=2;}
    if(fabs(a[1][2])>mx){mx=fabs(a[1][2]);p=1;q=2;}
    if(mx<1e-12) break;
    double phi=0.5*atan2(2*a[p][q], a[q][q]-a[p][p]), c=cos(phi),s=sin(phi);
    for(int k=0;k<3;k++){double kp=a[k][p],kq=a[k][q]; a[k][p]=c*kp-s*kq; a[k][q]=s*kp+c*kq;}
    for(int k=0;k<3;k++){double pk=a[p][k],qk=a[q][k]; a[p][k]=c*pk-s*qk; a[q][k]=s*pk+c*qk;}
    for(int k=0;k<3;k++){double kp=v[k][p],kq=v[k][q]; v[k][p]=c*kp-s*kq; v[k][q]=s*kp+c*kq;}
  }
  double e[3]={a[0][0],a[1][1],a[2][2]}; int o[3]={0,1,2};
  for(int i=0;i<3;i++)for(int j=i+1;j<3;j++) if(fabs(e[o[j]])<fabs(e[o[i]])){int t=o[i];o[i]=o[j];o[j]=t;}
  for(int i=0;i<3;i++){ ev[i]=e[o[i]]; for(int k=0;k<3;k++) V[k][i]=v[k][o[i]]; }
}

int hessian_sheet_detect(const f32 *vol, int nz, int ny, int nx,
                         const hess_params *p, f32 *sheetness, f32 *normal){
  size_t nynx=(size_t)ny*nx, n=(size_t)nz*nynx;
  static const f32 defs[3]={1.0f,2.0f,3.0f};
  const f32 *sig = (p&&p->sigmas)? p->sigmas : defs;
  int ns = (p&&p->sigmas)? p->nsigma : 3;
  float alpha = p? p->alpha : 0.5f;
  if (sheetness) for(size_t i=0;i<n;i++) sheetness[i]=0;

  f32 *sm=malloc(n*sizeof(f32));
  for(int si=0; si<ns; si++){
    double sigma=sig[si];
    gblur(vol,sm,nz,ny,nx,sigma);
    double g2=sigma*sigma;             // gamma-normalization (gamma=1)
    double c = (p && p->beta_c>0)? p->beta_c : 15.0;   // structure-strength cutoff
    #pragma omp parallel for schedule(static)
    for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
      int xm=CL(x-1,nx),xp=CL(x+1,nx),ym=CL(y-1,ny),yp=CL(y+1,ny),zm=CL(z-1,nz),zp=CL(z+1,nz);
      double cc=sm[IDX(z,y,x)];
      double Hxx=g2*(sm[IDX(z,y,xp)]-2*cc+sm[IDX(z,y,xm)]);
      double Hyy=g2*(sm[IDX(z,yp,x)]-2*cc+sm[IDX(z,ym,x)]);
      double Hzz=g2*(sm[IDX(zp,y,x)]-2*cc+sm[IDX(zm,y,x)]);
      double Hxy=g2*0.25*(sm[IDX(z,yp,xp)]-sm[IDX(z,yp,xm)]-sm[IDX(z,ym,xp)]+sm[IDX(z,ym,xm)]);
      double Hxz=g2*0.25*(sm[IDX(zp,y,xp)]-sm[IDX(zp,y,xm)]-sm[IDX(zm,y,xp)]+sm[IDX(zm,y,xm)]);
      double Hyz=g2*0.25*(sm[IDX(zp,yp,x)]-sm[IDX(zp,ym,x)]-sm[IDX(zm,yp,x)]+sm[IDX(zm,ym,x)]);
      double A[3][3]={{Hxx,Hxy,Hxz},{Hxy,Hyy,Hyz},{Hxz,Hyz,Hzz}};
      double ev[3],V[3][3]; eig3m(A,ev,V);
      double l1=ev[0],l2=ev[1],l3=ev[2];            // |l1|<=|l2|<=|l3|
      double S=sqrt(l1*l1+l2*l2+l3*l3);
      double resp=0.0;
      if (l3<0.0 && fabs(l3)>1e-9){                  // bright plate
        double Ra=fabs(l2)/fabs(l3);
        double plate=exp(-Ra*Ra/(2.0*(double)alpha*alpha));
        double strength=1.0-exp(-S*S/(2.0*c*c));
        resp=plate*strength;
      }
      size_t i=IDX(z,y,x);
      if (sheetness && resp>sheetness[i]){
        sheetness[i]=(f32)resp;
        if (normal){ // eigenvector of l3 (across-sheet)
          double nx_=V[0][2],ny_=V[1][2],nz_=V[2][2];
          double L=sqrt(nx_*nx_+ny_*ny_+nz_*nz_); if(L>1e-9){nx_/=L;ny_/=L;nz_/=L;}
          normal[3*i+0]=(f32)nx_; normal[3*i+1]=(f32)ny_; normal[3*i+2]=(f32)nz_;
        }
      }
    }
  }
  free(sm);
  return 0;
}
