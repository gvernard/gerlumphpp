#include "profile.hpp"


Profile::Profile(double pixSizePhys,double Rein,double width,double height){
  this->pixSizePhys = pixSizePhys;

  this->Nx = (int) ceil(width/pixSizePhys);
  this->Ny = (int) ceil(height/pixSizePhys);

  // Making sure that the profile dimensions in pixels are always even
  if( this->Nx%2 != 0 ){
    this->Nx += 1;
  }
  if( this->Ny%2 != 0 ){
    this->Ny += 1;
  }

  this->width  = Nx*pixSizePhys/Rein; // in Rein
  this->height = Ny*pixSizePhys/Rein; // in Rein

  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
}


void Profile::readProfile(const std::string filename){
  // read the profile from some file format, most likely fits, and interpolate on the appropriate resolution
  this->normalize();
}


void Profile::createGaussian(double gauss_width,double gauss_height,double incl,double orient){
  // still need to implement inclination and orientation

  double x  = 0.0;
  double y  = 0.0;
  double sx = 2*pow(gauss_width,2);
  double sy = 2*pow(gauss_height,2);
  
  for(int i=0;i<this->Ny/2;i++){
    y = (i - this->Ny/2)*this->pixSizePhys + this->pixSizePhys/2.0;
    for(int j=0;j<this->Nx/2;j++){
      x = (j - this->Nx/2)*this->pixSizePhys + this->pixSizePhys/2.0;
      
      this->data[i*this->Nx+j] = exp(-x*x/sx-y*y/sy);
    }
  }

  this->normalize();
}


void Profile::normalize(){
  double sum = 0.;
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      sum += this->data[i*Nx+j];
      //      std::cout << array[j*Nx+j] << std::endl;
    }
  }
  //  std::cout << sum << std::endl;
  
  for(int i=0;i<Ny;i++){
    for(int j=0;j<Nx;j++){
      this->data[i*Nx+j] /= sum;
    }
  }
}



