/* wind_poisson — INDEPENDENT winding solve for cross-validating sheet_sep3d (review P1).
 *
 * Solves the winding with a COMPLETELY DIFFERENT mechanism: the unwrap library's
 * winding_field_solve (a Poisson/Laplace diffusion warm-started from a contour-warped Archimedean
 * init), vs sheet_sep3d's radial-step-edge graph IRLS. Same archive, same region, same umbilicus,
 * same pitch -> if the two AGREE, the winding is method-independent (robust); if they diverge, the
 * disagreement localizes where one method is wrong. Dumps a _vol.f32 (same layout as sheet_sep3d's:
 * int32 {dz,dy,dx,z0,y0,x0} then dz*dy*dx float32, NaN off-material) so they compare directly.
 *
 * NOTE: winding_field produces a SPIRAL coordinate (r/pitch + theta/2pi); sheet_sep3d's is ~RADIAL
 * (constant along a sheet). To compare per-voxel, convert sheet_sep3d's w -> w + theta/2pi first
 * (scripts/compare_winding.py does this). The WRAP COUNT (span) is directly comparable.
 *
 *   wind_poisson ARCHIVE OUTBASE lod z0 y0 x0 dz dy dx pitch [iters=120] [anchor=0.02]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"
#include "annotate/umbilicus.h"
#include "unwrap/winding_field.h"
#include "segmentation/sheet_tensor.h"

static int air_threshold(const u8 *v, size_t n){
  long h[256]={0}; long tot=0; double allsum=0; long allnz=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x){ allsum+=x; allnz++; } if(x>=1&&x<=254){ h[x]++; tot++; } }
  if(tot<256){ int t=allnz?(int)(0.5*allsum/allnz+0.5):1; if(t<1)t=1; if(t>254)t=254; return t; }
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sumB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){ wB+=h[t]; if(wB==0)continue; double wF=(double)tot-wB; if(wF<=0)break;
    sumB+=(double)t*h[t]; double mB=sumB/wB,mF=(sum-sumB)/wF,btw=wB*wF*(mB-mF)*(mB-mF);
    if(btw>best){best=btw;thr=t;} }
  return thr;
}

int main(int argc,char**argv){
  if(argc<11){ fprintf(stderr,"usage: %s ARC OUTBASE lod z0 y0 x0 dz dy dx pitch [iters=120] [anchor=0.02]\n",argv[0]); return 2; }
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]);
  int dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]); double pitch=atof(argv[10]);
  int iters=argc>11?atoi(argv[11]):120; double anchor=argc>12?atof(argv[12]):0.02;
  int usenormals=argc>13?atoi(argv[13]):0;   // P2: forcing = div(sheet-normal/pitch) -> follow ACTUAL sheets

  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  // global umbilicus from coarse LOD5 (same as sheet_sep3d)
  int cl=5; double cs=(double)(1<<cl); int ccz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,ccz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse fail\n");return 1;}
  size_t ccn=(size_t)ccz*ccy*ccx; u8*ccm=malloc(ccn); for(size_t i=0;i<ccn;i++)ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,ccz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;} free(ccm);
  double ls=(double)(1<<lod);
  // subvol-local umbilicus: sample the global axis at NC control z, convert LOD5->subvol-local
  const int NC=9; umbilicus ul; ul.n=NC; ul.z=malloc(NC*sizeof(f32)); ul.y=malloc(NC*sizeof(f32)); ul.x=malloc(NC*sizeof(f32));
  for(int i=0;i<NC;i++){ double zz=(double)i*(dz-1)/(NC-1); double coarse_z=(z0+zz)*ls/cs;
    if(coarse_z<0)coarse_z=0; if(coarse_z>ccz-1)coarse_z=ccz-1; f32 ucy,ucx; umbilicus_center(&umb,(f32)coarse_z,&ucy,&ucx);
    ul.z[i]=(f32)zz; ul.y[i]=(f32)(ucy*cs/ls - y0); ul.x[i]=(f32)(ucx*cs/ls - x0); }
  fprintf(stderr,"subvol z%ld y%ld x%ld + %dx%dx%d; center z0=(%.0f,%.0f) zmid=(%.0f,%.0f) pitch=%.1f\n",
          z0,y0,x0,dz,dy,dx,ul.y[0],ul.x[0],ul.y[NC/2],ul.x[NC/2],pitch);

  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,dz,dy,dx); if(!v){fprintf(stderr,"subvol read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr?1:0; nmat+=mask[p]; }
  fprintf(stderr,"air<%d, %ld/%zu material\n",athr,nmat,nn);

  // P2 NORMAL FORCING: instead of the uniform-pitch radial assumption, drive the Poisson solve to
  // match the ACTUAL sheet geometry. Desired winding gradient g = n_out/pitch (n_out = sheet normal
  // oriented outward from the umbilicus); forcing = div(g) turns the Laplacian into a Poisson solve
  // whose level sets follow the real sheets, capturing variable density a uniform-pitch field misses.
  f32*forcing=NULL;
  if(usenormals){
    f32*nrm=malloc(3*nn*sizeof(f32)); f32*vf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++)vf[p]=v[p];
    st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=(f32)(pitch*0.12); if(sp.sigma_tensor<1)sp.sigma_tensor=1; if(sp.sigma_tensor>3)sp.sigma_tensor=3;
    st_sheet_detect(vf,dz,dy,dx,&sp,NULL,nrm); free(vf);
    f32*gx=calloc(nn,sizeof(f32)),*gy=calloc(nn,sizeof(f32)),*gz=calloc(nn,sizeof(f32));
    for(int z=0;z<dz;z++){ f32 cyz,cxz; umbilicus_center(&ul,(f32)z,&cyz,&cxz);
      for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x; if(!mask[p])continue;
        double nx=nrm[3*p+0],ny=nrm[3*p+1],nz=nrm[3*p+2];
        double rx=x-cxz,ry=y-cyz; if(nx*rx+ny*ry<0){nx=-nx;ny=-ny;nz=-nz;}   // orient outward
        gx[p]=(f32)(nx/pitch); gy[p]=(f32)(ny/pitch); gz[p]=(f32)(nz/pitch); } }
    forcing=calloc(nn,sizeof(f32));
    for(int z=1;z<dz-1;z++)for(int y=1;y<dy-1;y++)for(int x=1;x<dx-1;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x; if(!mask[p])continue;
      double div=0.5*((double)gx[p+1]-gx[p-1]) + 0.5*((double)gy[p+dx]-gy[p-dx]) + 0.5*((double)gz[p+(size_t)dy*dx]-gz[p-(size_t)dy*dx]);
      forcing[p]=(f32)div; }
    free(nrm);free(gx);free(gy);free(gz);
    fprintf(stderr,"normal forcing ON (div of sheet-normal/pitch)\n");
  }
  f32*wind=malloc(nn*sizeof(f32));
  winding_contour_warp(mask,dz,dy,dx,&ul,pitch,wind);   // Archimedean envelope init
  wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)pitch; wp.iters=iters; wp.omega=0.3f;
  wp.warm_start=1; wp.anchor_lambda=(f32)anchor; wp.forcing=forcing;   // anchor preserves spiral monodromy
  if(winding_field_solve(mask,dz,dy,dx,&ul,&wp,NULL,NULL,wind)){fprintf(stderr,"solve fail\n");return 1;}

  double wmin=1e30,wmax=-1e30; for(size_t p=0;p<nn;p++) if(mask[p]){ if(wind[p]<wmin)wmin=wind[p]; if(wind[p]>wmax)wmax=wind[p]; }
  fprintf(stderr,"POISSON winding: %.1f..%.1f (%.0f wraps) over %dx%dx%d\n",wmin,wmax,wmax-wmin,dz,dy,dx);

  char fn[700]; snprintf(fn,sizeof fn,"%s_vol.f32",base); FILE*f=fopen(fn,"wb");
  if(f){ int hdr[6]={dz,dy,dx,(int)z0,(int)y0,(int)x0}; fwrite(hdr,sizeof(int),6,f);
    float*wo=malloc(nn*sizeof(float)); for(size_t p=0;p<nn;p++) wo[p]=mask[p]?wind[p]:NAN;
    fwrite(wo,sizeof(float),nn,f); free(wo); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  umbilicus_free(&umb); free(ul.z);free(ul.y);free(ul.x); free(v);free(mask);free(wind); if(forcing)free(forcing); mca_close(arc);
  return 0;
}
