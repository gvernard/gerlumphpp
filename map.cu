#include "map.hpp"


Map::Map(std::string id){
  this->id = id;
  this->convolved = false;
  std::string file;


  // Read map metadata
  file = this->path + this->id + "/mapmeta.dat";
  std::ifstream myfile(file.c_str());
  myfile >> this->avgmu >> this->avgN;
  myfile >> this->Nx;
  myfile >> this->width;
  myfile >> this->k >> this->g >> this->s;
  myfile.close();
  this->Ny = this->Nx;
  this->height = this->width;
  this->pixSizeRein = this->width/this->Nx;


  // Read map data
  file = this->path + this->id + "/map.bin";

  FILE* ptr_myfile = fopen(file.data(),"rb");
  int* imap = (int*) calloc(this->Nx*this->Ny,sizeof(int));
  fread(imap,sizeof(int),this->Nx*this->Ny,ptr_myfile);
  fclose(ptr_myfile);
  
  //int (4 bytes) and cufftDoubleReal (8 bytes) do not have the same size, so there has to be a type cast
  double factor = fabs(this->avgmu/this->avgN);
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = (double) imap[i]*factor;
  }
  free(imap);
}


Map::Map(const Map& other){
  this->id     = other.id;
  this->k      = other.k;
  this->g      = other.g;
  this->s      = other.s;
  this->Nx     = other.Nx;
  this->Ny     = other.Ny;
  this->width  = other.width;
  this->height = other.height;
  this->avgmu  = other.avgmu;
  this->avgN   = other.avgN;
  this->pixSizeRein = other.pixSizeRein;
  this->pixSizePhys = other.pixSizePhys; // in units of [10^14 cm]
  this->convolved   = other.convolved;

  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = other.data[i];
  }
}


void Map::convolve(Profile* profile){
  cufftDoubleReal dum1,dum2;
  
  // Check if "kernel", which is a "profile" variable has the same dimension as the map

  //Fourier transform map
  cufftDoubleComplex* Fmap = (cufftDoubleComplex*) calloc(this->Nx*(this->Ny/2+1),sizeof(cufftDoubleComplex));
  myfft2d_r2c(this->Nx,this->Ny,this->data,Fmap);
  //Fourier transform kernel
  cufftDoubleComplex* Fkernel = (cufftDoubleComplex*) calloc(this->Nx*(this->Ny/2+1),sizeof(cufftDoubleComplex));
  myfft2d_r2c(this->Nx,this->Ny,profile->kernel,Fkernel);
  //Multiply kernel and map
  for(long i=0;i<this->Nx*(this->Ny/2+1);i++){
    dum1 = (cufftDoubleReal) (Fmap[i].x*Fkernel[i].x - Fmap[i].y*Fkernel[i].y);
    dum2 = (cufftDoubleReal) (Fmap[i].x*Fkernel[i].y + Fmap[i].y*Fkernel[i].x);
    Fmap[i].x = dum1;
    Fmap[i].y = dum2;
  }
  //Inverse Fourier transform
  //  cmap = (cufftDoubleReal*) calloc(this->res*this->res,sizeof(cufftDoubleReal));
  //  myfft2d_c2r(this->res,this->res,Fmap,cmap);
  myfft2d_c2r(this->Nx,this->Ny,Fmap,this->data);

  free(Fmap);
  free(Fkernel);


  //Normalize convolved map
  double norm = (double) this->Nx*this->Ny;
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->data[i*this->Nx+j] /= norm;
    }
  }

  this->convolved = true;
}


Mpd* Map::getFullMpd(){
  if( this->convolved ){
    std::cout << "Map is convolved. Has to be in ray counts." << std::endl;
    return NULL;
  } else {
    thrust::device_vector<double> bins;
    thrust::device_vector<int> counts;
    thrust::device_vector<double> data(this->data,this->data+this->Nx*this->Ny);
    thrust::sort(data.begin(),data.end());
    int num_bins = thrust::inner_product(data.begin(),data.end()-1,data.begin()+1,int(1),thrust::plus<int>(),thrust::not_equal_to<double>());
    counts.resize(num_bins);
    bins.resize(num_bins);
    thrust::reduce_by_key(data.begin(),data.end(),thrust::constant_iterator<int>(1),counts.begin(),bins.begin());
    thrust::host_vector<int> hcounts(counts);
    thrust::host_vector<double> hbins(bins);
    
    Mpd* theMpd = new Mpd(hcounts.size());
    for(unsigned int i=0;i<hcounts.size();i++){
      theMpd->counts[i] = (double) hcounts[i]/(this->Nx*this->Ny);
      theMpd->bins[i]   = (double) hbins[i]*this->avgmu/this->avgN;
    }
    return theMpd; 
  }
}


Mpd* Map::getBinnedMpd(int Nbins){
  double min  = log10(0.02);
  double max  = log10(200);
  double dbin = (max-min)/Nbins;
  double* bins = (double*) calloc(Nbins,sizeof(double));
  for(int i=0;i<Nbins;i++){
    bins[i] = pow(10,min+(i+1)*dbin);
  }

  thrust::device_vector<int> counts(Nbins);
  thrust::device_vector<double> dbins(bins,bins+Nbins);
  thrust::device_vector<double> data(this->data,this->data+this->Nx*this->Ny);
  thrust::sort(data.begin(), data.end());
  thrust::upper_bound(data.begin(),data.end(),dbins.begin(),dbins.end(),counts.begin());
  thrust::adjacent_difference(counts.begin(),counts.end(),counts.begin());
  thrust::host_vector<int> hcounts(counts);

  Mpd* theMpd = new Mpd(hcounts.size());
  for(unsigned int i=0;i<hcounts.size();i++){
    theMpd->counts[i] = (double) hcounts[i]/(this->Nx*this->Ny);
    theMpd->bins[i]   = (double) bins[i]*this->avgmu/this->avgN;
  }
  return theMpd;
}


int Map::myfft2d_r2c(int Nx,int Ny,cufftDoubleReal* data,cufftDoubleComplex* Fdata){
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


int Map::myfft2d_c2r(int Nx, int Ny, cufftDoubleComplex* Fdata, cufftDoubleReal* data){
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

