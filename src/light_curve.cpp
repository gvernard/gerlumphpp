#include "light_curve.hpp"

LightCurveCollection::LightCurveCollection(int Ncurves,EffectiveMap* emap){
  this->Ncurves = Ncurves;
  this->A = (point*) calloc(Ncurves,sizeof(point));
  this->B = (point*) calloc(Ncurves,sizeof(point));
  this->lightCurves = (LightCurve*) calloc(Ncurves,sizeof(LightCurve));
  this->pixSizePhys = emap->pixSizePhys;
  this->Nx = emap->Nx;
  this->Ny = emap->Ny;
  this->emap = emap;
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
  this->pixSizePhys = other.pixSizePhys;
  this->Nx = other.Nx;
  this->Ny = other.Ny;
  this->emap = other.emap;
}

void LightCurveCollection::setEmap(EffectiveMap* emap){
  this->emap = emap;
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

void LightCurveCollection::createVelocityLocations(int seed,double tmax,std::vector<double> v,std::vector<double> phi){
  srand48(seed);
  double d2r = 0.017453; // degrees to radians
  point A;
  point B;
  double len,lenx,leny,phi_rad;

  for(int i=0;i<this->Ncurves;i++){
    len = tmax*v[i];
    phi_rad = phi[i]*d2r;
    lenx = len*cos(phi_rad);
    leny = len*sin(phi_rad);
    
    while(true){
      A.x = drand48()*this->Nx;
      A.y = drand48()*this->Ny;
    
      B.x = A.x + lenx;
      B.y = A.y + leny;
      
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

std::vector<int> LightCurveCollection::checkLengthFull(){
  // all comparisons are made in pixels and the index of problematic light curves is returned
  std::vector<int> indices;
  for(int i=0;i<this->Ncurves;i++){
    int Dj = floor(this->B[i].x - this->A[i].x);
    int Di = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    if( Lmax < this->full_sampling_lower_limit ){
      indices.push_back(i);
    }
  }
  return indices;
}

std::vector<int> LightCurveCollection::checkLength(double v,double tmax){
  // all comparisons are made in pixels and the index of problematic light curves is returned
  std::vector<int> indices;
  int lmax = floor(v*tmax*this->factor/this->pixSizePhys);
  for(int i=0;i<this->Ncurves;i++){
    int Dj = floor(this->B[i].x - this->A[i].x);
    int Di = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    if( lmax >= Lmax ){
      indices.push_back(i);
    }
  }
  return indices;
}

std::vector<int> LightCurveCollection::checkSampling(double v,double dt){
  std::vector<int> indices;
  double dl = v*dt*this->factor/this->pixSizePhys;
  for(int i=0;i<this->Ncurves;i++){
    int Dj = floor(this->B[i].x - this->A[i].x);
    int Di = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    if( Lmax/dl < this->sampling_lower_limit ){
      indices.push_back(i);
    }
  }
  return indices;
}

std::vector<int> LightCurveCollection::checkSampling(double v,std::vector<double> t){
  std::vector<int> indices;
  for(int i=0;i<this->Ncurves;i++){
    int Dj = floor(this->B[i].x - this->A[i].x);
    int Di = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));

    int j=0;
    for(j=0;j<t.size();j++){
      double lj = v*t[j]*this->factor/this->pixSizePhys;
      if( lj > Lmax ){
	break;
      }
    }
    if( j < this->sampling_lower_limit ){
      indices.push_back(i);
    }
  }
  return indices;
}



void LightCurveCollection::extractFull(){
  this->type = "full";
  for(int i=0;i<this->Ncurves;i++){
    int Dj   = floor(this->B[i].x - this->A[i].x);
    int Di   = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    int dl   = 1;
    int Nsamples = Lmax;
    double phi = atan2(Di,Dj);

    this->lightCurves[i].Nsamples = Nsamples;
    this->lightCurves[i].t  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].m  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].dm = (double*) calloc(Nsamples,sizeof(double));

    std::vector<double> length(Nsamples); // in pixels (but can be decimal number as well)
    for(int k=0;k<Nsamples;k++){
      length[k] = k;
      this->lightCurves[i].t[k] = length[k]*this->pixSizePhys; // in [10^14 cm]
    }

    sampleLightCurve(i,length,phi);
  }
}

void LightCurveCollection::extractSampled(double v,double dt,double tmax){
  std::vector<double> vvec(this->Ncurves);
  std::fill(vvec.begin(),vvec.end(),v); 
  extractSampled(vvec,dt,tmax);
}

void LightCurveCollection::extractSampled(std::vector<double> v,double dt,double tmax){
  this->type = "sampled";
  for(int i=0;i<this->Ncurves;i++){
    int Dj   = floor(this->B[i].x - this->A[i].x);
    int Di   = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    int dl   = v[i]*dt*this->factor/this->pixSizePhys;
    int Nsamples = floor(Lmax/dl);
    double phi = atan2(Di,Dj);

    this->lightCurves[i].Nsamples = Nsamples;
    this->lightCurves[i].t  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].m  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].dm = (double*) calloc(Nsamples,sizeof(double));

    std::vector<double> length(Nsamples); // in pixels (but can be decimal number as well)
    for(int k=0;k<Nsamples;k++){
      length[k] = k*dl;
      this->lightCurves[i].t[k] = length[k]*this->pixSizePhys; // in [10^14 cm]
    }

    sampleLightCurve(i,length,phi);
  }
}

void LightCurveCollection::extractStrategy(double v,std::vector<double> t){
  std::vector<double> vvec(this->Ncurves);
  std::fill(vvec.begin(),vvec.end(),v); 
  extractStrategy(vvec,t);
}

void LightCurveCollection::extractStrategy(std::vector<double> v,std::vector<double> t){
  this->type = "strategy";
  for(int i=0;i<this->Ncurves;i++){
    int Dj   = floor(this->B[i].x - this->A[i].x);
    int Di   = floor(this->B[i].y - this->A[i].y);
    int Lmax = floor(hypot(Di,Dj));
    double phi = atan2(Di,Dj);

    std::vector<double> length; 
    for(int j=0;j<t.size();j++){
      double len = v[i]*t[j]*this->factor/this->pixSizePhys;
      if( len <= Lmax ){
	length.push_back(len);
      } else {
	break;
      }
    }
    int Nsamples = length.size();

    this->lightCurves[i].Nsamples = Nsamples;
    this->lightCurves[i].t  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].m  = (double*) calloc(Nsamples,sizeof(double));
    this->lightCurves[i].dm = (double*) calloc(Nsamples,sizeof(double));

    for(int k=0;k<Nsamples;k++){
      this->lightCurves[i].t[k] = length[k]*this->pixSizePhys; // in [10^14 cm]
    }

    sampleLightCurve(i,length,phi);
  }
}

void LightCurveCollection::sampleLightCurve(int index,std::vector<double> length,double phi){
  double m,dm;
  int Nsamples = length.size();

  for(int k=0;k<Nsamples;k++){
    double xk = this->A[index].x + length[k]*cos(phi);
    double yk = this->A[index].y + length[k]*sin(phi);
  
    interpolatePlane(xk,yk,m,dm);
    
    this->lightCurves[index].m[k]  = m;
    this->lightCurves[index].dm[k] = dm;
  }
}

void LightCurveCollection::interpolatePlane(double xk,double yk,double& m,double& dm){
  int j0 = floor(xk);
  int i0 = floor(yk);

  double dx  = (double) fabs(j0 - xk);
  double dy  = (double) fabs(i0 - yk);
  double ddx = (double) (1.0 - dx);
  double ddy = (double) (1.0 - dy);
  double w00 = dx*dy;
  double w10 = dy*ddx;
  double w01 = dx*ddy;
  double w11 = ddx*ddy;

  double f00 = this->emap->data[i0*this->emap->Nx+j0];
  double f10 = this->emap->data[i0*this->emap->Nx+j0+1];
  double f01 = this->emap->data[(i0+1)*this->emap->Nx+j0];
  double f11 = this->emap->data[(i0+1)*this->emap->Nx+j0+1];
  m = f00*w00 + f10*w10 + f01*w01 + f11*w11;

  double df00 = sqrt(f00*fabs(this->emap->avgmu/this->emap->avgN));
  double df10 = sqrt(f10*fabs(this->emap->avgmu/this->emap->avgN));
  double df01 = sqrt(f01*fabs(this->emap->avgmu/this->emap->avgN));
  double df11 = sqrt(f11*fabs(this->emap->avgmu/this->emap->avgN));
  dm = df00*w00 + df10*w10 + df01*w01 + df11*w11;
}


void LightCurveCollection::writeCurves(const std::string prefix){
  for(int i=0;i<this->Ncurves;i++){
    std::string fname = prefix + std::to_string(i) + ".dat";
    this->lightCurves[i].writeData(fname);
  }
}


void LightCurveCollection::writeCurvesDegraded(const std::string prefix,const std::string degraded){
  if( this->type == "full" || this->type == "sampled" ){
    try {
      if( degraded == "byte" || degraded == "int16" ){
	for(int i=0;i<this->Ncurves;i++){
	  std::string fname = prefix + std::to_string(i) + "_b.bin";
	  this->lightCurves[i].writeDegraded(fname,degraded); // "degraded" has to be "byte" or "int16"
	}
      } else {
	throw "Only allowed options for a degraded output with a full or regularly sampled light curve are: 'byte' and 'int16'";
      }
    } catch(const char* msg){
      std::cout << msg << std::endl;
    }
  } else if( this->type == "strategy" ){
    try {
      if( degraded == "bytebyte" ){
	for(int i=0;i<this->Ncurves;i++){
	  std::string fname = prefix + std::to_string(i) + "_bb.bin";
	  this->lightCurves[i].writeDegraded(fname,"byte","byte");
	}
      } else if( degraded == "int16byte" ){
	for(int i=0;i<this->Ncurves;i++){
	  std::string fname = prefix + std::to_string(i) + "_ib.bin";
	  this->lightCurves[i].writeDegraded(fname,"int16","byte");
	}
      } else if( degraded == "int16int16" ){
	for(int i=0;i<this->Ncurves;i++){
	  std::string fname = prefix + std::to_string(i) + "_ii.bin";
	  this->lightCurves[i].writeDegraded(fname,"int16","int16");
	}
      } else {
	throw "Only allowed options for a degraded output with an irregularly sampled light curve are: 'bytebyte', 'int16byte', and 'int16int16'";
      }
    } catch(const char* msg){
      std::cout << msg << std::endl;
    }
  }
}

void LightCurve::writeData(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Nsamples;i++){
    fprintf(fh,"%11.6e %11.6e %11.6e\n",this->t[i],this->m[i],this->dm[i]);
  }
  fclose(fh);
}

void LightCurve::writeDegraded(const std::string filename,std::string m_type){

}

void LightCurve::writeDegraded(const std::string filename,std::string t_type,std::string m_type){

}
