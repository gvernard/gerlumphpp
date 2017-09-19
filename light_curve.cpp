#include "light_curve.hpp"


LightCurveCollection::LightCurveCollection(int Ncurves,EffectiveMap* emap){
  this->Ncurves = Ncurves;
  this->A = (point*) calloc(Ncurves,sizeof(point));
  this->B = (point*) calloc(Ncurves,sizeof(point));
  this->lightCurves = (LightCurve*) calloc(Ncurves,sizeof(LightCurve));
  this->pixSizePhys = emap->pixSizePhys;
  this->Nx = emap->Nx;
  this->Ny = emap->Ny;
}

LightCurveCollection::LightCurveCollection(const LightCurveCollection& other){
  this->Ncurves = other.Ncurves;
  this->A = (point*) calloc(other.Ncurves,sizeof(point));
  this->B = (point*) calloc(other.Ncurves,sizeof(point));

  for(int i=0;i<other.Ncurves;i++){
    this->A[i] = other.A[i];
    this->B[i] = other.B[i];
  }
  
  this->lightCurves = (LightCurve*) calloc(other.Ncurves,sizeof(LightCurve));
  this->pixSizePhys = pixSizePhys;
}

void LightCurveCollection::createRandomLocations(int seed,int maxLen){
  srand48(seed);
  point A;
  point B;
  double ranx,rany,len;

  for(int i=0;i<this->Ncurves;i++){
    A.x = drand48()*this->Nx;
    A.y = drand48()*this->Ny;
    
    while(true){
      // point in a random direction within the map
      ranx = drand48()*this->Nx;
      rany = drand48()*this->Ny;
      
      // set the actual length and end point for the light curve
      len = hypot((A.x-ranx),(A.y-rany));
      B.x = A.x + (ranx-A.x)*maxLen/len;
      B.y = A.y + (rany-A.y)*maxLen/len;
      
      //Is the new end point still within the (effective) map?
      if ( (B.x < 0) || (B.x >= this->Nx) || (B.y < 0) || (B.y >= this->Ny) ){
	//	printf("%5d%5d problem\n",(int)round(loc.xo),(int)round(loc.yo));
	continue;
      } else {
	break;
      }
    }
    
    this->A[i] = A;
    this->B[i] = B;
  }
}

void LightCurveCollection::writeLocations(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Ncurves;i++){
    fprintf(fh,"%5d%5d%5d%5d\n",(int)round(this->A[i].x),(int)round(this->A[i].y),(int)round(this->B[i].x),(int)round(this->B[i].y));
  }
  fclose(fh);
}

