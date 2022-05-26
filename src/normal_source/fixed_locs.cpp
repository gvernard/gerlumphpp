#include "fixed_locs.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>



FixedLocationCollection::FixedLocationCollection(int Nlocs,EffectiveMap* emap){
  this->Nlocs = Nlocs;
  this->Nx = emap->Nx;
  this->Ny = emap->Ny;
  this->emap = emap;
  this->A  = (point*) calloc(Nlocs,sizeof(point));
  this->m  = (double*) calloc(Nlocs,sizeof(double));
  this->dm = (double*) calloc(Nlocs,sizeof(double));
}

FixedLocationCollection::FixedLocationCollection(int Nlocs,int Nx,int Ny){
  // Need to call setEmap manually if this constructor is used
  this->Nlocs = Nlocs;
  this->Nx = Nx;
  this->Ny = Ny;
  this->A  = (point*) calloc(Nlocs,sizeof(point));
  this->m  = (double*) calloc(Nlocs,sizeof(double));
  this->dm = (double*) calloc(Nlocs,sizeof(double));
}

FixedLocationCollection::FixedLocationCollection(const FixedLocationCollection& other){
  this->Nlocs = other.Nlocs;
  this->type  = other.type;
  this->A  = (point*) calloc(other.Nlocs,sizeof(point));
  this->m  = (double*) calloc(other.Nlocs,sizeof(double));
  this->dm = (double*) calloc(other.Nlocs,sizeof(double));

  for(int i=0;i<other.Nlocs;i++){
    this->A[i] = other.A[i];
  }

  this->Nx = other.Nx;
  this->Ny = other.Ny;
  this->emap = other.emap;
}

void FixedLocationCollection::setEmap(EffectiveMap* emap){
  this->emap = emap;
}

void FixedLocationCollection::createRandomLocations(int seed){
  this->type = "random";
  srand48(seed);
  point A;
  
  for(int i=0;i<this->Nlocs;i++){
    A.x = drand48()*this->Nx;
    A.y = drand48()*this->Ny;
    this->A[i] = A;
  }
}

void FixedLocationCollection::createGridLocations(){
  this->type = "grid";
  try {
    double test = sqrt(this->Nlocs);
    if( test != floor(test) ){
      throw "Need to reset the number of locations to a square number (representing a square grid). Locations uninitialized.";
    }
    int len = floor(test);
    point A;
    int dpix_x = this->Nx/len;
    int dpix_y = this->Nx/len;
    int offset_x = dpix_x/2.0;
    int offset_y = dpix_y/2.0;
    for(int i=0;i<len;i++){
      A.y = floor(offset_y + i*dpix_y);      
      for(int j=0;j<len;j++){
	A.x = floor(offset_x + j*dpix_x);
	this->A[i*len+j] = A;
      }
    }
  } catch(const char* msg) {
    std::cout << msg << " " << std::endl;
  }
}

double FixedLocationCollection::checkOverlap(int profRadius){
  double pi = 3.14159;
  point P1,P2;
  if( this->type == "grid" ){
    P1 = this->A[0];
    P2 = this->A[1];
    double d = hypot(P2.x-P1.x,P2.y-P1.y);
    if( d < profRadius ){
      double dum = d/(2.0*profRadius);
      double theta = asin( sqrt((1.0+dum)*(1.0-dum)) )/2.0;
      double overlap = (theta-sin(theta))/pi;
      return overlap;
    } else {
      return 0.0;
    }
  } else if( this->type == "random" ){
    double avg_overlap = 0;
    for( int i=0;i<this->Nlocs-1;i++){
      P1 = this->A[i];
      for( int j=i+1;j<this->Nlocs;j++){
	P2 = this->A[j];

	double d = hypot(P2.x-P1.x,P2.y-P1.y);
	if( d < profRadius ){
	  double dum = d/(2.0*profRadius);
	  double theta = asin( sqrt((1.0+dum)*(1.0-dum)) )/2.0;
	  avg_overlap += (theta-sin(theta));
	}
      }
    }
    int N_unique = (this->Nlocs-1)*this->Nlocs/2;
    avg_overlap /= N_unique*pi;
    return avg_overlap;
  } else {
    return 0.0;
  }
}


void FixedLocationCollection::writeLocations(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Nlocs;i++){
    fprintf(fh,"%5d%5d\n",(int)round(this->A[i].x),(int)round(this->A[i].y));
  }
  fclose(fh);
}

void FixedLocationCollection::extract(){
  for(int k=0;k<this->Nlocs;k++){
    int j = floor(this->A[k].x);
    int i = floor(this->A[k].y);
    this->m[k]  = this->emap->data[i*this->Nx+j];
    this->dm[k] = sqrt(this->m[k]*fabs(this->emap->avgmu/this->emap->avgN));
  }  
}

void FixedLocationCollection::writeData(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Nlocs;i++){
    fprintf(fh,"%11.6e %11.6e\n",this->m[i],this->dm[i]);
  }
  fclose(fh);
}

