#include "kernel.hpp"

Kernel::Kernel(int map_Nx,int map_Ny,Profile* profile){
  int fNx = profile->Nx;   // full profile width in pixels
  int fNy = profile->Ny;   // full profile height in pixels
  int hNx = profile->Nx/2; // half profile width in pixels
  int hNy = profile->Ny/2; // half profile height in pixels

  this->data = (double*) calloc(map_Nx*map_Ny,sizeof(double));
  for(int i=0;i<hNy;i++){
    for(int j=0;j<hNx;j++){
      this->data[i*map_Nx+j]                                = profile->data[hNy*fNx+hNx+i*fNx+j];
      this->data[map_Nx-hNx+i*map_Nx+j]                     = profile->data[hNy*fNx+i*fNx+j];
      this->data[map_Nx*(map_Ny-hNy)+i*map_Nx+j]            = profile->data[hNx+fNx*i+j];
      this->data[map_Nx*(map_Ny-hNy)+map_Nx-hNx+i*map_Nx+j] = profile->data[fNx*i+j];
    }
  }

  this->Nx  = map_Nx;
  this->Ny  = map_Ny;
  this->hNx = profile->Nx/2;
  this->hNy = profile->Ny/2;
}

