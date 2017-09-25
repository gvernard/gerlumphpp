#include "fixed_locs.hpp"

FixedLocationCollection::FixedLocationCollection(int Nlocs,EffectiveMap* emap){
  this->Nlocs = Nlocs;
  this->Nx = emap->Nx;
  this->Ny = emap->Ny;
  this->emap = emap;
  this->A  = (point*) calloc(Nlocs,sizeof(point));
  this->m  = (double*) calloc(Nlocs,sizeof(double));
  this->dm = (double*) calloc(Nlocs,sizeof(double));
}

FixedLocationCollection::FixedLocationCollection(const FixedLocationCollection& other){
  this->Nlocs = other.Nlocs;
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
  srand48(seed);
  point A;
  
  for(int i=0;i<this->Nlocs;i++){
    A.x = drand48()*this->Nx;
    A.y = drand48()*this->Ny;
    this->A[i] = A;
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

