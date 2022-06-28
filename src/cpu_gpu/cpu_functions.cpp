#include "magnification_map.hpp"

#include <fftw3.h>
#include <algorithm>
#include <cstring>

using namespace gerlumph;

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
  Mpd theMpd(0);
  try {
    if( this->convolved ){
      throw "Map is convolved. It has to be in ray counts in order to use this function.";
    }

    long int N = this->Nx*this->Ny;
    double* tmp = (double*) malloc(N*sizeof(double));
    std::memcpy(tmp,this->data,N*sizeof(double));
    std::sort(tmp,tmp+N);

    int i = 0;
    std::vector<double> bins;
    std::vector<int> counts;
    while( i < N ){
      double val0 = tmp[i];
      int j = 0;
      while( i < N ){
	i++;
	if( tmp[i] == val0 ){
	  j++;
	} else {
	  break;
	}
      }
      bins.push_back(val0);
      counts.push_back(j);
    }
    free(tmp);

    theMpd.reset(bins.size());
    for(int i=0;i<bins.size();i++){
      theMpd.counts[i] = (double) (counts[i])/(double) (N);
      theMpd.bins[i] = bins[i];
    }
    
  } catch(const char* msg){
    std::cout << msg << std::endl;
  }
  return theMpd;
}

Mpd MagnificationMap::getBinnedMpd(int Nbins){
  // creating bins which are evenly spaced in log space
  double min = 0.02;
  double max = 200;

  double logmin  = log10(min);
  double logmax  = log10(max);
  double logdbin = (logmax-logmin)/Nbins;
  double* bins   = (double*) calloc(Nbins,sizeof(double));
  for(int i=0;i<Nbins;i++){
    bins[i] = pow(10,logmin+(i+1)*logdbin);
  }

 
  long int N = this->Nx*this->Ny;
  double* tmp = (double*) malloc(N*sizeof(double));
  std::memcpy(tmp,this->data,N*sizeof(double));
  std::sort(tmp,tmp+N); // sorts in ascending order

  
  std::vector<int> counts(Nbins);
  std::fill(counts.begin(),counts.end(),0);
  int i = 0;
  while( tmp[i] < min ){
    i++;
  }
  int j = 0;
  while( i < N && j < Nbins ){
    if( tmp[i] <= bins[j] ){
      counts[j] += 1;
      i++;
    } else {
      j++;
    }
  }
  free(tmp);

  Mpd theMpd(Nbins);
  for(int i=0;i<Nbins;i++){
    if( counts[i] == 0 ){
      theMpd.counts[i] = 1.e-8;
    } else {
      theMpd.counts[i] = (double) (counts[i])/(double) (N);
    }
    theMpd.bins[i] = bins[i];
  }  
  free(bins);
  
  return theMpd;
}
