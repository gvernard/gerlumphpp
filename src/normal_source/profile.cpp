#include "profile.hpp"
#include "rectGrid.hpp"

#include <cmath>
#include <iostream>
#include <string>

#include <CCfits/CCfits>

using namespace gerlumph;

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
double BaseProfile::getSize(std::map<std::string,std::string> pars,double lrest){
  // If the temperature profile is parameteric or SS then calculate the half-light radius.
  // If it is simply a vector of given half-light radii, then just lrest is the rhalf anyway.
  // If it is a custom profile, then calculate the half-light radius numerically.
  double rhalf;
  if( pars["type"] == "parametric" ){
    double r0 = std::stof(pars["r0"]);
    double l0 = std::stof(pars["l0"]);
    double nu = std::stof(pars["nu"]);
    rhalf = BaseProfile::sizeParametric(r0,l0,nu,lrest);
  } else if( pars["type"] == "ss_disc" ){
    double mbh  = std::stof(pars["mbh"]);
    double fedd = std::stof(pars["fedd"]);
    double eta  = std::stof(pars["eta"]);
    rhalf = BaseProfile::sizeSS(mbh,fedd,eta,lrest);
  } else if( pars["type"] == "vector" ){
    rhalf = lrest;
  } else if( pars["type"] == "custom" ){
    rhalf = 1.0;
  } else {
    rhalf = 0.0; // throw an exception
  }
  return rhalf;
}
double BaseProfile::sizeParametric(double r0,double l0,double nu,double lrest){
  // r0 in [10^14 cm], l0 and lrest in [nm]
  double r = r0*pow(lrest/l0,nu);
  return r; // in [10^14 cm]
}
double BaseProfile::sizeSS(double mbh,double fedd,double eta,double lrest){
  double a = pow(lrest,4.0);
  double b = pow(mbh,2.0);
  double c = fedd/eta;
  double r = 0.0097*pow(a*b*c,1.0/3.0); // in [10^14 cm]
  return r;
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
double UniformDisc::setByHalfRadius(double rhalf){
  return sqrt(2.0)*rhalf;
}
double UniformDisc::getHalfRadius(){
  return this->R/sqrt(2.0);
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
double Gaussian::setByHalfRadius(double rhalf){
  return rhalf/1.18;
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
double GaussianHole::setByHalfRadius(double rhalf,double Rin){
  double dum = (pow(rhalf,2.0) - pow(Rin,2.0))/(2.0*log(2.0));
  return sqrt(dum);
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
double ThermalHole::setByHalfRadius(double rhalf){
  return rhalf/4.0;
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
double Wavy::setByHalfRadius(double rhalf,int node){
  return node*M_PI/(2.0*rhalf);
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
double Exponential::setByHalfRadius(double rhalf){
  return 0.4106*rhalf;
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
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  int Nxx = (int) image.axis(0);
  int Nyy = (int) image.axis(1);

  double xmax = Nxx*profPixSizePhys;
  double ymax = Nyy*profPixSizePhys;
  RectGrid input_profile(Nxx,Nyy,0,xmax,0,ymax,filename);
  
  int newNx = (int) floor(xmax/pixSizePhys);
  int newNy = (int) floor(ymax/pixSizePhys);
  RectGrid new_profile = input_profile.embeddedNewGrid(newNx,newNy);

  this->Nx = new_profile.Nx;
  this->Ny = new_profile.Ny;
  this->width = xmax;
  this->height = ymax;
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->data[i*this->Nx+j] = new_profile.z[i*this->Nx+j];
    }
  }
  
  normalize();
}

double Custom::getHalfRadius(){
  double center_x = this->width/2.0;
  double center_y = this->height/2.0;

  // Calculate the distance of each pixel from the center.
  std::vector<double> radii(Nx*Ny);
  std::vector<double> values(Nx*Ny);
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      radii[i*this->Nx+j] = hypot(i*this->pixSizePhys-center_y,j*this->pixSizePhys-center_x);
      values[i*this->Nx+j] = this->data[i*this->Nx+j];
    }
  }
  // Sort both distances and values by the distance.
  auto p = sort_permutation(radii);
  radii = apply_permutation(radii,p);
  values = apply_permutation(values,p);

  // Now we have a list of all the pixels ordered by distance from the center.
  // Just go through the list and add the elements.
  // The half-light radius will be at the index where the sum is equal to 0.5 (the profile is normalized in the constructor).
  double sum = 0.0;
  double rhalf = 3.3;
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      sum += values[i*Nx+j];
      if( sum > 0.5 ){
	rhalf = radii[i*Nx+j];
	break;
      }
    }
  }
  
  return rhalf;
}

std::vector<std::size_t> Custom::sort_permutation(const std::vector<double>& vec){
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),[&](std::size_t i, std::size_t j){ return vec[i]<vec[j]; });
  return p;
}

std::vector<double> Custom::apply_permutation(const std::vector<double>& vec,const std::vector<std::size_t>& p){
  std::vector<double> sorted_vec(vec.size());
  std::transform(p.begin(), p.end(), sorted_vec.begin(),[&](std::size_t i){ return vec[i]; });
  return sorted_vec;
}
