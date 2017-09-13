class Profile {
public:
  double size;   // [in cm]
  double incl;   // inclination angle in degrees
  double orient; // orientation angle in degrees
  double* kernel;


  void normalizeBinned(int Nx,int Ny,double* array){
    double sum = 0.;
    for(int j=0;j<Ny;j++){
      for(int i=0;i<Nx;i++){
	sum += array[j*Nx+i];
	//      std::cout << array[j*Nx+j] << std::endl;
      }
    }
    //  std::cout << sum << std::endl;
    
    for(int j=0;j<Ny;j++){
      for(int i=0;i<Nx;i++){
	array[j*Nx+i] /= sum;
      }
    }
  }

  void setKernel(int bNx,int bNy,double* binned,int Nx,int Ny,cufftDoubleReal* kernel){
    for(long k=0;k<Nx*Ny;k++){
      kernel[k] = 0;
    }
    
    for(int j=0;j<bNy;j++){
      for(int i=0;i<bNx;i++){
	kernel[j*Nx+i]                    = (cufftDoubleReal) binned[bNy*2*bNx+bNx+j*2*bNx+i];
	kernel[Nx-bNx+j*Nx+i]             = (cufftDoubleReal) binned[bNy*2*bNx+j*2*bNx+i];
	kernel[Nx*(Ny-bNy)+j*Nx+i]        = (cufftDoubleReal) binned[bNx+2*bNx*j+i];
	kernel[Nx*(Ny-bNy)+Nx-bNx+j*Nx+i] = (cufftDoubleReal) binned[2*bNx*j+i];
      }
    }
  }

  
};
