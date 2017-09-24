#include "profile.hpp"


Profile::Profile(double pixSizePhys,double Rein,double width,double height){
  this->imageType = "profile";
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


Profile::Profile(double pixSizePhys,const std::string filename,double profPixSizePhys){
  this->imageType = "profile";
  this->pixSizePhys = pixSizePhys;

  
  // Read the profile from .fits format
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  int Nxx = image.axis(0);
  int Nyy = image.axis(1);
  double* input = (double*) calloc(Nxx*Nyy,sizeof(double));

  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int i=0;i<Nyy;i++){
    for(int j=0;j<Nxx;j++){
      input[i*Nxx+j] = (double) tmp[(Nyy-i-1)*Nxx+j];
    }
  }

  // Depending on the profile and map pixel sizes in physical units, select whether to interpolate or bin the profile to match the map
  if( this->pixSizePhys <= profPixSizePhys ){
    interpolateProfile(Nxx,Nyy,input,profPixSizePhys);
  } else {
    binProfile(Nxx,Nyy,input,profPixSizePhys);
  }
  free(input);

  this->normalize();
}


void Profile::interpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys){
  // Decide on the profile width and height in pixels based on the input profile
  double dresx = (Nxx-1)*profPixSizePhys/this->pixSizePhys;
  this->Nx = floor(dresx);
  double xoffset = (dresx - this->Nx)/2.0;
  if( this->Nx%2 != 0 ){
    this->Nx -= 1;
    xoffset += 0.5;
  }
  double dresy = (Nyy-1)*profPixSizePhys/this->pixSizePhys;
  this->Ny = floor(dresy);
  double yoffset = (dresy - this->Ny)/2.0;
  if( this->Ny%2 != 0 ){
    this->Ny -= 1;
    yoffset += 0.5;
  }
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));

  double x,y,xp,yp,dx,dy,ddx,ddy,w00,w10,w01,w11,f00,f10,f01,f11;
  int ii,jj;
  for(int i=0;i<this->Ny;i++){
    y  = yoffset+i*this->pixSizePhys;
    ii = floor( y/profPixSizePhys );
    yp = ii*profPixSizePhys;
    dy = y - yp;
    ddy = (1.0 - dy);

    for(int j=0;j<this->Nx;j++){
      x  = xoffset+j*this->pixSizePhys;
      jj = floor( x/profPixSizePhys );
      xp = jj*profPixSizePhys;
      dx = x - xp;
      ddx = (1.0 - dx);

      w00 = dx*dy;
      w10 = dy*ddx;
      w01 = dx*ddy;
      w11 = ddx*ddy;

      f00 = input[ii*Nxx+jj];
      f10 = input[ii*Nxx+jj+1];
      f01 = input[(ii+1)*Nxx+jj];
      f11 = input[(ii+1)*Nxx+jj+1];

      this->data[i*this->Nx+j] = f00*w00 + f10*w10 + f01*w01 + f11*w11;
    }
  }
}

void Profile::binProfile(int Nxx,int Nyy,double* input,double profPixSizePhys){
  double dresx = Nxx*profPixSizePhys/this->pixSizePhys;
  this->Nx = ceil(dresx);
  double xoffset = (this->Nx - dresx)/2.0;
  if( this->Nx%2 != 0 ){
    this->Nx += 1;
    xoffset += 0.5;
  }
  double dresy = Nyy*profPixSizePhys/this->pixSizePhys;
  this->Ny = ceil(dresy);
  double yoffset = (this->Ny - dresy)/2.0;
  if( this->Ny%2 != 0 ){
    this->Ny += 1;
    yoffset += 0.5;
  }
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));

  int* bin_counts = (int*) calloc(this->Nx*this->Ny,sizeof(int));
  int ii,jj;
  double x,y;
  for(int i=0;i<Nyy;i++){
    y  = yoffset+i*profPixSizePhys;
    ii = floor( y/this->pixSizePhys );

    for(int j=0;j<Nxx;j++){
      x  = xoffset+j*profPixSizePhys;
      jj = floor( x/this->pixSizePhys );

      this->data[ii*this->Nx+jj] += input[i*Nxx+j];
      bin_counts[ii*this->Nx+jj]++;
    }
  }

  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = this->data[i]/((double) bin_counts[i]);
  }
  free(bin_counts);
}






void Profile::createGaussian(double gauss_width,double gauss_height,double incl,double orient){
  // still need to implement inclination and orientation
  double x  = 0.0;
  double y  = 0.0;
  double sx = 2*pow(gauss_width,2);
  double sy = 2*pow(gauss_height,2);
  
  for(int i=0;i<this->Ny;i++){
    y = (i - this->Ny/2)*this->pixSizePhys + this->pixSizePhys/2.0;
    for(int j=0;j<this->Nx;j++){
      x = (j - this->Nx/2)*this->pixSizePhys + this->pixSizePhys/2.0;
      this->data[i*this->Nx+j] = exp(-x*x/sx-y*y/sy);
    }
  }

  this->normalize();
}


void Profile::createUniDisc(double radius,double incl,double orient){
  // still need to implement inclination and orientation
  double x,y,r;
  for(int i=0;i<this->Ny;i++){
    y = (i - this->Ny/2)*this->pixSizePhys + this->pixSizePhys/2.0;
    for(int j=0;j<this->Nx;j++){
      x = (j - this->Nx/2)*this->pixSizePhys + this->pixSizePhys/2.0;
      r = hypot(x,y);
      if( r < radius ){
	this->data[i*this->Nx+j] = 1.0;
      } else {
	this->data[i*this->Nx+j] = 0.0;
      }
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
      //      std::cout << this->data[i*Nx+j] << std::endl;
   }
  }
}



/*
void Profile::writeProfileFits(const std::string filename){
  long naxis    = 2;
  long naxes[2] = {(long) this->Nx,(long) this->Ny};
  long Ntot = (long) this->Nx*this->Ny;
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::auto_ptr<CCfits::FITS> pFits(0);
  pFits.reset( new CCfits::FITS("!"+filename,FLOAT_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) this->Ny);
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  //Need Ni and Nj as index counters to flip image
  std::valarray<float> array(Ntot);
  long count = 0;
  for(int j=0;j<this->Nx;j++){
    for(int i=0;i<this->Ny;i++){
      array[(this->Nx-1-j)*this->Ny+i] = (float) (this->data[count]);
      count++;
    }
  }
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,array);
  //  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
  //  std::cout << pFits->pHDU() << std::endl;
}
*/
