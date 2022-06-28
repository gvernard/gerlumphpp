#include "fitsInterface.hpp"

#include <CCfits/CCfits>
#include <valarray>

using namespace gerlumph;

void FitsInterface::readFits(int Nx,int Ny,double* z,const std::string filepath){
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  
  int NNy = image.axis(0);
  int NNx = image.axis(1);
  if( NNy != Ny || NNx != Nx ){
    // throw an exception
  }
  
  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int i=0;i<Ny;i++){
    for(int j=0;j<Nx;j++){
      //z[i*Nx+j] = tmp[(Ny-i-1)*Nx+j];
      z[i*Nx+j] = tmp[i*Nx+j];
    }
  }
}

void FitsInterface::writeFits(int Nx,int Ny,double* z,const std::string filepath){
  std::vector<std::string> dum;
  FitsInterface::writeFits(Nx,Ny,z,dum,dum,dum,filepath);
}

void FitsInterface::writeFits(int Nx,int Ny,double* z,std::vector<std::string> key,std::vector<std::string> value,std::vector<std::string> description,const std::string filepath){
  long naxis    = 2;
  long naxes[2] = {(long) Nx,(long) Ny};
  long Ntot = (long) Nx*Ny;
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::unique_ptr<CCfits::FITS> pFits( new CCfits::FITS("!"+filepath,FLOAT_IMG,naxis,naxes) );

  std::vector<long> extAx(2,(long) Ny);
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  std::valarray<float> array(Ntot);
  long count = 0;
  for(int i=0;i<Ny;i++){
    for(int j=0;j<Nx;j++){
      //array[(Ni-1-i)*Nj+j] = (float) z[count];
      array[i*Nx+j] = (float) z[count];
      count++;
    }
  }
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,array);
  //  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 

  for(int k=0;k<key.size();k++){
    pFits->pHDU().addKey(key[k],value[k],description[k]);
  }
  pFits->pHDU().write(fpixel,Ntot,array); 
  //  std::cout << pFits->pHDU() << std::endl;
}
