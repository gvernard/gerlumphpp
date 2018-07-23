#include "magnification_map.hpp"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <thrust/inner_product.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <cufft.h>

// These two prototype functions are used only to compile this source code.
// The 'static' keyword causes their visibility to be limited to the translation unit (this single .cu or .cpp file where they are defined).
static int myfft2d_r2c(int Nx, int Ny, cufftDoubleReal* data, cufftDoubleComplex* Fdata);
static int myfft2d_c2r(int Nx, int Ny, cufftDoubleComplex* Fdata, cufftDoubleReal* data);


void MagnificationMap::convolve(Kernel* kernel,EffectiveMap* emap){
  cufftDoubleReal dum1,dum2;
  
  // Check if "kernel", which is a "profile" variable has the same dimension as the map
  
  //Fourier transform map
  cufftDoubleComplex* Fmap = (cufftDoubleComplex*) calloc(this->Nx*(this->Ny/2+1),sizeof(cufftDoubleComplex));
  myfft2d_r2c(this->Nx,this->Ny,this->data,Fmap);
  //Fourier transform kernel
  cufftDoubleComplex* Fkernel = (cufftDoubleComplex*) calloc(this->Nx*(this->Ny/2+1),sizeof(cufftDoubleComplex));
  myfft2d_r2c(this->Nx,this->Ny,kernel->data,Fkernel);
  //Multiply kernel and map
  for(long i=0;i<this->Nx*(this->Ny/2+1);i++){
    dum1 = (cufftDoubleReal) (Fmap[i].x*Fkernel[i].x - Fmap[i].y*Fkernel[i].y);
    dum2 = (cufftDoubleReal) (Fmap[i].x*Fkernel[i].y + Fmap[i].y*Fkernel[i].x);
    Fmap[i].x = dum1;
    Fmap[i].y = dum2;
  }
  //Inverse Fourier transform
  cufftDoubleReal* cmap = (cufftDoubleReal*) calloc(this->Nx*this->Ny,sizeof(cufftDoubleReal));
  myfft2d_c2r(this->Nx,this->Ny,Fmap,cmap);

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
    thrust::device_vector<int> counts;
    thrust::device_vector<double> bins;
    thrust::device_vector<double> data(this->data,this->data+this->Nx*this->Ny);
    thrust::sort(data.begin(),data.end());

    int num_bins = thrust::inner_product(data.begin(),data.end()-1,data.begin()+1,int(1),thrust::plus<int>(),thrust::not_equal_to<double>());
    counts.resize(num_bins);
    bins.resize(num_bins);

    thrust::reduce_by_key(data.begin(),data.end(),thrust::constant_iterator<int>(1),bins.begin(),counts.begin());
    thrust::host_vector<int> hcounts(counts);
    thrust::host_vector<double> hbins(bins);
    
    theMpd.reset(num_bins);
    for(unsigned int i=0;i<hcounts.size();i++){
      theMpd.counts[i] = (double) (hcounts[i])/(double) (this->Nx*this->Ny);
      theMpd.bins[i]   = ((double) (hbins[i]));
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

  thrust::device_vector<int>    counts(Nbins);
  thrust::device_vector<double> dbins(bins,bins+Nbins);
  thrust::device_vector<double> data(this->data,this->data+this->Nx*this->Ny);
  thrust::sort(data.begin(),data.end());

  // For the following lines to work I need to compile using the flag: --expt-extended-lambda
  //  auto getLog10LambdaFunctor = [=]  __device__ (double x) {return log10(x);};
  //  thrust::transform(data.begin(),data.end(),data.begin(),getLog10LambdaFunctor);

  double range[2] = {min,max};
  thrust::device_vector<double> drange(range,range+2);
  thrust::device_vector<int>    dirange(2);
  thrust::lower_bound(data.begin(),data.end(),drange.begin(),drange.end(),dirange.begin());
  thrust::host_vector<int> hirange(dirange);
  //  std::cout << hirange[0] << " " << hirange[1] << std::endl;

  thrust::upper_bound(data.begin() + hirange[0],data.begin() + hirange[1],dbins.begin(),dbins.end(),counts.begin());
  //  thrust::upper_bound(data.begin(),data.end(),dbins.begin(),dbins.end(),counts.begin());
  thrust::adjacent_difference(counts.begin(),counts.end(),counts.begin());
  thrust::host_vector<int>    hcounts(counts);

  Mpd theMpd(hcounts.size());
  for(unsigned int i=0;i<hcounts.size();i++){
    theMpd.counts[i] = (double) (hcounts[i]) /(double) (this->Nx*this->Ny);
    theMpd.bins[i]   = (double) bins[i];
  }
  free(bins);
  return theMpd;
}


int myfft2d_r2c(int Nx,int Ny,cufftDoubleReal* data,cufftDoubleComplex* Fdata){
  cufftHandle plan;
  cufftDoubleReal* data_GPU;
  cufftDoubleComplex* Fdata_GPU;

  //allocate and transfer memory to the GPU
  cudaMalloc( (void**) &data_GPU, Nx*Ny*sizeof(cufftDoubleReal));
  cudaMemcpy( data_GPU, data, Nx*Ny*sizeof(cufftDoubleReal), cudaMemcpyHostToDevice);
  cudaMalloc( (void**) &Fdata_GPU, Nx*(Ny/2+1)*sizeof(cufftDoubleComplex));

  //do the fourier transform on the GPU
  cufftPlan2d(&plan,Nx,Ny,CUFFT_D2Z);
  cufftExecD2Z(plan, data_GPU, Fdata_GPU);
  //  cudaDeviceSynchronize();
  cudaThreadSynchronize();
  cufftDestroy(plan);

  //copy back results
  cudaMemcpy(Fdata, Fdata_GPU, Nx*(Ny/2+1)*sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
  cudaFree(data_GPU);
  cudaFree(Fdata_GPU);

  return 0;
}


int myfft2d_c2r(int Nx, int Ny, cufftDoubleComplex* Fdata, cufftDoubleReal* data){
  cufftHandle plan;
  cufftDoubleComplex* Fdata_GPU;
  cufftDoubleReal* data_GPU;
  
  //allocate and transfer memory
  cudaMalloc((void**) &Fdata_GPU, Nx*(Ny/2+1)*sizeof(cufftDoubleComplex));
  cudaMemcpy(Fdata_GPU, Fdata, Nx*(Ny/2+1)*sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
  cudaMalloc((void**) &data_GPU, Nx*Ny*sizeof(cufftDoubleReal));

  //do the inverse fourier transform on the GPU
  cufftPlan2d(&plan,Nx,Ny,CUFFT_Z2D) ;
  cufftExecZ2D(plan, Fdata_GPU, data_GPU);
  //  cudaDeviceSynchronize();
  cudaThreadSynchronize();
  cufftDestroy(plan);
  
  //copy back results
  cudaMemcpy(data, data_GPU, Nx*Ny*sizeof(cufftDoubleReal), cudaMemcpyDeviceToHost);
  cudaFree(data_GPU);
  cudaFree(Fdata_GPU);
  
  return 0;
}
