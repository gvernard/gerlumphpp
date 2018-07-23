#include "magnification_map.hpp"

#include <fftw3.h>

void MagnificationMap::convolve(Kernel* kernel,EffectiveMap* emap){
  double dum1,dum2;
  fftw_plan p1;
  
  //Fourier transform map
  fftw_complex* Fmap    = (fftw_complex*) fftw_malloc(this->Nx*this->Ny*sizeof(fftw_complex));
  p1 = fftw_plan_dft_r2c_2d(this->Nx,this->Ny,this->data,Fmap,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  //Fourier transform kernel
  fftw_complex* Fkernel = (fftw_complex*) fftw_malloc(this->Nx*this->Ny*sizeof(fftw_complex));
  p1 = fftw_plan_dft_r2c_2d(this->Nx,this->Ny,kernel->data,Fkernel,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  //Multiply kernel and map
  for(int i=0;i<this->Nx*this->Ny;i++) {
    dum1 = Fmap[i][0]*Fkernel[i][0] - Fmap[i][1]*Fkernel[i][1];
    dum2 = Fmap[i][0]*Fkernel[i][1] + Fmap[i][1]*Fkernel[i][0];
    Fmap[i][0] = dum1;
    Fmap[i][1] = dum2;
  }
  //Inverse Fourier transform
  double* cmap = (double*)  malloc(this->Nx*this->Ny*sizeof(double));
  p1 = fftw_plan_dft_c2r_2d(this->Nx,this->Ny,Fmap,cmap,FFTW_ESTIMATE);
  fftw_execute(p1);
  fftw_destroy_plan(p1);
  
  free(Fmap);
  free(Fkernel);
  
  //Normalize convolved map and crop to emap
  double norm = (double) (this->Nx*this->Ny);
  for(int i=0;i<emap->Ny;i++){
    for(int j=0;j<emap->Nx;j++){
      emap->data[i*emap->Nx+j] = (double) (cmap[emap->top*this->Nx+emap->left+i*this->Nx+j]/norm);
    }
  }
  
  this->convolved = true;
  free(cmap);
}

Mpd MagnificationMap::getFullMpd(){
}

Mpd MagnificationMap::getBinnedMpd(int Nbins){
}
