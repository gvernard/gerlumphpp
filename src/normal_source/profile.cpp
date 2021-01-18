#include "profile.hpp"

#include <cmath>
#include <iostream>
#include <string>

#include <CCfits/CCfits>



//////////////////////// CLASS IMPLEMENTATION: BaseProfile ////////////////////////
///////////////////////////////////////////////////////////////////////////////////
BaseProfile::BaseProfile(double pixSizePhys,double incl,double orient){
  this->imageType = "profile";
  this->pixSizePhys = pixSizePhys;
  this->incl = incl;
  this->orient = orient;
}
BaseProfile::BaseProfile(const BaseProfile& other){
  this->imageType = other.imageType;
  this->pixSizePhys = other.pixSizePhys;
  this->incl = other.incl;
  this->orient = other.orient;
  this->Nx = other.Nx;
  this->Ny = other.Ny;
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = other.data[i];
  }
}
double BaseProfile::sizeParametric(parsParametric pars,double lrest){
  // r0 in [10^14 cm], l0 and lrest in [nm]
  double r = pars.r0*pow(lrest/pars.l0,pars.nu);
  return r; // in [10^14 cm]
}
double BaseProfile::sizeSS(parsSSdisc pars,double lrest){
  double a = pow(lrest,4.0);
  double b = pow(pars.mbh,2.0);
  double c = pars.fedd/pars.eta;
  double r = 0.0097*pow(a*b*c,1.0/3.0); // in [10^14 cm]
}
void BaseProfile::project(int N,double* x,double* y){
  double fac = M_PI/180.0; // degrees to radians conversion
  double costheta = cos(fac*this->incl); // cos(incl)
  // orientation happens clockwise around the z axis, so I have to add pi/2 and change sign to have it from x-axis counter-clockwise
  double cosphi = cos(fac*(-this->orient)); // cos(orient) 
  double sinphi = sin(fac*(-this->orient)); // sin(orient)

  double newx,newy;
  for(int i=0;i<N;i++){
    //    newx = costheta*cosphi*x[i] - sinphi*y[i];
    //    newy = costheta*sinphi*x[i] + cosphi*y[i];
    newx = x[i]*cosphi/costheta + y[i]*sinphi/costheta;
    newy = -sinphi*x[i] + cosphi*y[i];
    x[i] = newx;
    y[i] = newy;
  }
}
void BaseProfile::createGrid(double* x,double* y){
  // Nx and Ny must be set beforehand
  // first create the orthogonal x,y grid of the image
  double dumx,dumy;
  for(int i=0;i<this->Ny;i++){
    dumy = (i - this->Ny/2)*this->pixSizePhys + this->pixSizePhys/2.0;
    for(int j=0;j<this->Nx;j++){
      dumx = (j - this->Nx/2)*this->pixSizePhys + this->pixSizePhys/2.0;
      x[i*this->Nx+j] = dumx;
      y[i*this->Nx+j] = dumy;
    }
  }
  // then rotate the x,y grid according to the inclination and orientation angles
  this->project(this->Nx*this->Ny,x,y);
}
void BaseProfile::normalize(){
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
void BaseProfile::makeEven(int& N){
  if( N%2 != 0 ){
    N += 1;
  }
}



//////////////////////// CLASS IMPLEMENTATION: UniformDisc ////////////////////////
///////////////////////////////////////////////////////////////////////////////////
UniformDisc::UniformDisc(double pixSizePhys,double R,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->R = R;
  generateValues();
  normalize();
}
UniformDisc::UniformDisc(double pixSizePhys,parsParametric pars,double lrest,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->R = sizeParametric(pars,lrest);
  generateValues();
  normalize();
}
UniformDisc::UniformDisc(double pixSizePhys,parsSSdisc pars,double lrest,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->R = sizeSS(pars,lrest);
  generateValues();
  normalize();
}
double UniformDisc::getHalfRadius(){
  return 0.707*this->R;
}
void UniformDisc::generateValues(){
  this->Nx = 2 * (3+(int) ceil(this->R/this->pixSizePhys)); // x2 the disc radius + 3 pixels
  this->Ny = 2 * (3+(int) ceil(this->R/this->pixSizePhys)); // x2 the disc radius + 3 pixels
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data   = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    if( r < this->R ){
      this->data[i] = 1.0;
    } else {
      this->data[i] = 0.0;
    }
  }

  free(x);
  free(y);
}



//////////////////////// CLASS IMPLEMENTATION: Gaussian ////////////////////////
////////////////////////////////////////////////////////////////////////////////
Gaussian::Gaussian(double pixSizePhys,double sdev,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sdev;
  generateValues();
  normalize();
}
Gaussian::Gaussian(double pixSizePhys,parsParametric pars,double lrest,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sizeParametric(pars,lrest);
  generateValues();
  normalize();
}
Gaussian::Gaussian(double pixSizePhys,parsSSdisc pars,double lrest,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sizeSS(pars,lrest);
  generateValues();
  normalize();
}
double Gaussian::getHalfRadius(){
  return 1.18*this->sdev;
}
void Gaussian::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // width is equal to x4 the half light radius
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // height is equal to x4 the half light radius
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data  = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  double s = 2*pow(this->sdev,2);
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    this->data[i] = exp(-r*r/s);
  }

  free(x);
  free(y);
}



//////////////////////// CLASS IMPLEMENTATION: GaussianHole //////////////////
//////////////////////////////////////////////////////////////////////////////
GaussianHole::GaussianHole(double pixSizePhys,double sdev,double Rin,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sdev;
  this->Rin = Rin;
  generateValues();
  normalize();
}
GaussianHole::GaussianHole(double pixSizePhys,parsParametric pars,double lrest,double Rin,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sizeParametric(pars,lrest);
  this->Rin = Rin;
  generateValues();
  normalize();
}
GaussianHole::GaussianHole(double pixSizePhys,parsSSdisc pars,double lrest,double Rin,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sizeSS(pars,lrest);
  this->Rin = Rin;
  generateValues();
  normalize();
}
double GaussianHole::getHalfRadius(){
  double dum = 2.0*log(2.0)*pow(this->sdev,2.0) + pow(this->Rin,2.0);
  return sqrt(dum);
}
void GaussianHole::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // width is equal to x4 the half light radius
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // height is equal to x4 the half light radius
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data  = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  double s = 2*pow(this->sdev,2);
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    if( r < this->Rin ){
      this->data[i] = 0.0;
    } else {
      this->data[i] = exp(-r*r/s);
    }
  }

  free(x);
  free(y);
}

//////////////////////// CLASS IMPLEMENTATION: ThermalHole /////////////////////
////////////////////////////////////////////////////////////////////////////////
ThermalHole::ThermalHole(double pixSizePhys,double Rin,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->Rin = Rin;
  generateValues();
  normalize();
}
double ThermalHole::getHalfRadius(){
  return 4.0*this->Rin;
}
void ThermalHole::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // width is equal to x4 the half light radius
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // height is equal to x4 the half light radius
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data  = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  double fac = 3.0*this->Rin/(2.0*M_PI);
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    if( r > this->Rin ){
      this->data[i] = (1.0 - sqrt(this->Rin/r))*fac/pow(r,3.0);
    } else {
      this->data[i] = 0.0;
    }
  }

  free(x);
  free(y);
}

//////////////////////// CLASS IMPLEMENTATION: Wavy ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Wavy::Wavy(double pixSizePhys,double a,double b,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->a = a;
  this->b = b;
  this->node = 3; // the profile extends to the 3rd node of the sine wave.
  generateValues();
  normalize();
}
double Wavy::getHalfRadius(){
  return this->node*M_PI/(2.0*this->b);
}
void Wavy::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // width is equal to x4 the half light radius
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // height is equal to x4 the half light radius
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data  = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    if( r < 2*Rhalf ){
      this->data[i] = pow(sin(this->b*r),2)/(this->a*r);
    } else {
      this->data[i] = 0.0;
    }
  }

  free(x);
  free(y);
}

//////////////////////// CLASS IMPLEMENTATION: Exponential ////////////////////////
///////////////////////////////////////////////////////////////////////////////////
Exponential::Exponential(double pixSizePhys,double sigma,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sigma = sigma;
  generateValues();
  normalize();
}
double Exponential::getHalfRadius(){
  return this->sigma/0.4106;
}
void Exponential::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // x2 the disc radius + 3 pixels
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // x2 the disc radius + 3 pixels
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data   = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // Inner radius to truncate the profile because it is unbound for r->0, set to be the radius where the brightness is ten times the value at rhalf
  double v_inner = 10.0/(exp(pow(Rhalf/this->sigma,0.75)) - 1.0);
  double r_inner = this->sigma*pow(log(1.0+1.0/v_inner),4.0/3.0);

  // set the values
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    if( r < r_inner ){
      this->data[i] = v_inner;
    } if( r < 2*Rhalf ){
      this->data[i] = 1.0/(exp(pow(r/this->sigma,0.75)) - 1.0);
    } else {
      this->data[i] = 0.0;
    }
  }

  free(x);
  free(y);
}



//////////////////////// CLASS IMPLEMENTATION: Custom ////////////////////////
//////////////////////////////////////////////////////////////////////////////
Custom::Custom(double pixSizePhys,const std::string filename,double profPixSizePhys,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
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

  normalize();
}
void Custom::interpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys){
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
  this->data   = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

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
void Custom::binProfile(int Nxx,int Nyy,double* input,double profPixSizePhys){
  double dresx = (Nxx)*profPixSizePhys/this->pixSizePhys;
  this->Nx = ceil(dresx);
  double xoffset = fabs(this->Nx - dresx)*this->pixSizePhys/2.0;
  if( this->Nx%2 != 0 ){
    this->Nx += 1;
    xoffset += 0.5;
  }
  double dresy = (Nyy)*profPixSizePhys/this->pixSizePhys;
  this->Ny = ceil(dresy);
  double yoffset = fabs(this->Ny - dresy)*this->pixSizePhys/2.0;
  if( this->Ny%2 != 0 ){
    this->Ny += 1;
    yoffset += 0.5;
  }
  this->data   = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

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
      bin_counts[ii*this->Nx+jj] += 1;
    }
  }

  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = this->data[i]/((double) bin_counts[i]);
  }
  free(bin_counts);
}



//////////////////////// CLASS IMPLEMENTATION: Gaussian lamp-post ////////////////////////
////////////////////////////////////////////////////////////////////////////////
GaussianLP::GaussianLP(double pixSizePhys,double sdev,int t0,double* f,double incl,double orient) : BaseProfile(pixSizePhys,incl,orient) {
  this->sdev = sdev;
  this->t0 = t0;
  this->f = f;
  generateValues();
  normalize();
}
double GaussianLP::getHalfRadius(){
  return 1.18*this->sdev;
}
void GaussianLP::generateValues(){
  double Rhalf = this->getHalfRadius();
  this->Nx = (int) ceil(4.0*Rhalf/this->pixSizePhys); // width is equal to x4 the half light radius
  this->Ny = (int) ceil(4.0*Rhalf/this->pixSizePhys); // height is equal to x4 the half light radius
  makeEven(this->Nx);
  makeEven(this->Ny);
  this->data  = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width  = this->pixSizePhys*this->Nx;
  this->height = this->pixSizePhys*this->Ny;

  // create x,y grid
  double* x = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  double* y = (double*) malloc(this->Nx*this->Ny*sizeof(double));
  this->createGrid(x,y);

  // set the values
  double s = 2*pow(this->sdev,2);
  for(int i=0;i<this->Nx*this->Ny;i++){
    double r = hypot(x[i],y[i]);
    this->data[i] = exp(-r*r/s);
  }

  free(x);
  free(y);
}
