#include "light_curve.hpp"

#include <cstdlib>
#include <cmath>
#include <fstream>

LightCurveCollection::LightCurveCollection(int Ncurves){
  this->Ncurves = Ncurves;
  this->A = (point*) calloc(Ncurves,sizeof(point));
  this->B = (point*) calloc(Ncurves,sizeof(point));
  this->lightCurves = new LightCurve*[Ncurves];
}

LightCurveCollection::LightCurveCollection(int Ncurves,MagnificationMap* emap){
  this->Ncurves = Ncurves;
  this->A = (point*) calloc(Ncurves,sizeof(point));
  this->B = (point*) calloc(Ncurves,sizeof(point));
  this->lightCurves = new LightCurve*[Ncurves];
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
  
  this->lightCurves = new LightCurve*[Ncurves];
  for(int i=0;i<this->Ncurves;i++){
    this->lightCurves[i] = other.lightCurves[i];
  }
  this->pixSizePhys = other.pixSizePhys;
  this->Nx = other.Nx;
  this->Ny = other.Ny;
  this->emap = other.emap;
}

void LightCurveCollection::setEmap(MagnificationMap* emap){
  // If emap is reset to something with different resolution, I need to rescale A and B coordinates
  this->pixSizePhys = emap->pixSizePhys;
  this->Nx = emap->Nx;
  this->Ny = emap->Ny;
  this->emap = emap;
}

void LightCurveCollection::createRandomLocations(int seed,int maxLen){
  srand48(seed);
  point A;
  point B;
  double ranx,rany,phi;

  for(int i=0;i<this->Ncurves;i++){
    A.x = drand48()*(this->Nx - 1);
    A.y = drand48()*(this->Ny - 1);
    
    while(true){
      // point in a random direction within the map
      ranx = drand48()*(this->Nx - 1);
      rany = drand48()*(this->Ny - 1);
      
      // set the actual length and end point for the light curve
      phi = atan2(rany-A.y,ranx-A.x);
      B.x = A.x + cos(phi)*(maxLen);
      B.y = A.y + sin(phi)*(maxLen);
      
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

void LightCurveCollection::createOrientedRandomLocations(int seed,int maxLen,double angle){
  // "angle" is in degress measured counter-clockwise from the positive x-axis of the usual cartesian system.
  // The shear on the maps is oriented along this positive x-axis.
  srand48(seed);
  point A;
  point B;
  double fac = 0.017453293; // degrees to radians
  double angle_rad = angle*fac;
  
  for(int i=0;i<this->Ncurves;i++){
 
    while(true){
      // original starting point
      A.x = drand48()*(this->Nx - 1);
      A.y = drand48()*(this->Ny - 1);

      // point in the given direction within the map
      B.x = A.x + cos(angle_rad)*maxLen;
      B.y = A.y + sin(angle_rad)*maxLen;
      
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
  point A;
  point B;
  double len,lenx,leny,phi_rad;

  for(int i=0;i<this->Ncurves;i++){
    len = tmax*v[i]*this->vfactor/this->pixSizePhys; // has to be in pixel units
    phi_rad = phi[i]*this->d2r;
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

void LightCurveCollection::createVelocityLocations(int seed,double tmax,std::vector<double> v,std::vector<double> phi,double phig){
  // Same as above but it takes into account the rotation due to phig (external shear angle) to make sure the end point B is still inside the effective map.
  // Otherwise it changes starting point A and does it again.
  srand48(seed);
  point A;
  point B;
  double len,lenx,leny,phi_rad;

  for(int i=0;i<this->Ncurves;i++){
    len = tmax*v[i]*this->vfactor/this->pixSizePhys; // has to be in pixel units
    phi_rad = (phi[i] - phig - 90)*this->d2r; // the maps are at 90-degrees with respect to the shear orientation 2*phi = phig
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

void LightCurveCollection::createVelocityLocations(int seed,double tmax,std::vector<double> v,std::vector<double> phi,std::vector<double> phig){
  // Same as above but it takes into account a vector of rotations phig.
  // This function sets only starting point A.

  srand48(seed);
  point A;
  point B;
  double len,lenx,leny,phi_rad;

  for(int i=0;i<this->Ncurves;i++){
    len = tmax*v[i]*this->vfactor/this->pixSizePhys; // has to be in pixel units
    
    while(true){
      A.x = drand48()*this->Nx;
      A.y = drand48()*this->Ny;

      int false_counter = 0;
      for(int j=0;j<phig.size();j++){
	phi_rad = (phi[i] - phig[j] - 90)*this->d2r; // the maps are at 90-degrees with respect to the shear orientation 2*phi = phig
	lenx = len*cos(phi_rad);
	leny = len*sin(phi_rad);
	B.x = A.x + lenx;
	B.y = A.y + leny;
	//Is the new end point still within the (effective) map?
	if ( (B.x < 0) || (B.x >= this->Nx) || (B.y < 0) || (B.y >= this->Ny) ){
	  false_counter++;
	}
      }
      if( false_counter > 0 ){
	continue;
      } else {
	break;
      }
    }
    
    this->A[i] = A;
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
  int lmax = floor(v*tmax*this->vfactor/this->pixSizePhys);
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
  double dl = v*dt*this->vfactor/this->pixSizePhys;
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
      double lj = v*t[j]*this->vfactor/this->pixSizePhys;
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
  double dl  = this->pixSizePhys;
  for(int i=0;i<this->Ncurves;i++){
    double Dx   = this->B[i].x - this->A[i].x;
    double Dy   = this->B[i].y - this->A[i].y;
    int Nsamples = ceil(hypot(Dy,Dx)+0.5);
    double phi = atan2(Dy,Dx);

    this->lightCurves[i] = new LightCurve(Nsamples);

    std::vector<double> length(Nsamples); // in pixels (but can be decimal number as well)
    for(int k=0;k<Nsamples;k++){
      length[k] = k;
      this->lightCurves[i]->t[k] = length[k]*dl; // in [10^14 cm]
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
    int dl   = v[i]*dt*this->vfactor/this->pixSizePhys;
    int Nsamples = floor(Lmax/dl);
    double phi = atan2(Di,Dj);

    this->lightCurves[i] = new LightCurve(Nsamples);

    std::vector<double> length(Nsamples); // in pixels (but can be decimal number as well)
    for(int k=0;k<Nsamples;k++){
      length[k] = k*dl;
      this->lightCurves[i]->t[k] = length[k]*this->pixSizePhys; // in [10^14 cm]
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
    int Dj   = ceil(this->B[i].x - this->A[i].x);
    int Di   = ceil(this->B[i].y - this->A[i].y);
    int Lmax = ceil(hypot(Di,Dj));
    double phi = atan2(Di,Dj);

    std::vector<double> length; 
    for(int j=0;j<t.size();j++){
      double len = v[i]*t[j]*this->vfactor/this->pixSizePhys;
      if( len <= Lmax ){
	length.push_back(len);
      } else {
	break;
      }
    }
    int Nsamples = length.size();

    this->lightCurves[i] = new LightCurve(Nsamples);

    for(int k=0;k<Nsamples;k++){
      //this->lightCurves[i].t[k] = length[k]*this->pixSizePhys; // in [10^14 cm]
      this->lightCurves[i]->t[k] = t[k]; // in [days]
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
    
    this->lightCurves[index]->m[k]  = m;
    this->lightCurves[index]->dm[k] = dm;
  }
}

void LightCurveCollection::interpolatePlane(double xk,double yk,double& m,double& dm){
  int j0 = floor(xk);
  int i0 = floor(yk);

  double dx  = (double) fabs(j0 - xk);
  double dy  = (double) fabs(i0 - yk);
  double ddx = (double) (1.0 - dx);
  double ddy = (double) (1.0 - dy);
  double w00 = ddx*ddy;
  double w10 = ddy*dx;
  double w01 = ddx*dy;
  double w11 = dx*dy;

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


void LightCurveCollection::writeCurves(const std::string path,const std::string suffix){
  for(int i=0;i<this->Ncurves;i++){
    std::string fname = suffix + std::to_string(i) + ".dat";
    this->lightCurves[i]->writeData(path,fname);
  }
}
template<typename mType> void LightCurveCollection::writeCurvesDegraded(const std::string path,const std::string suffix){
  for(int i=0;i<this->Ncurves;i++){
    this->lightCurves[i]->writeDegraded<mType>(path,suffix + std::to_string(i));
  }
}
template<typename mType,typename tType> void LightCurveCollection::writeCurvesDegraded(const std::string path,const std::string suffix){
  for(int i=0;i<this->Ncurves;i++){
    this->lightCurves[i]->writeDegraded<mType,tType>(path,suffix + std::to_string(i));
  }
}
template<typename mType,typename tType,typename eType> void LightCurveCollection::writeCurvesDegraded(const std::string path,const std::string suffix){
  for(int i=0;i<this->Ncurves;i++){
    this->lightCurves[i]->writeDegraded<mType,tType,eType>(path,suffix + std::to_string(i));
  }
}

template void LightCurveCollection::writeCurvesDegraded<unsigned char>(const std::string path,const std::string suffix);
template void LightCurveCollection::writeCurvesDegraded<unsigned short int>(const std::string path,const std::string suffix);

template void LightCurveCollection::writeCurvesDegraded<unsigned char,unsigned char>(const std::string path,const std::string suffix);
template void LightCurveCollection::writeCurvesDegraded<unsigned char,unsigned short int>(const std::string path,const std::string suffix);

template void LightCurveCollection::writeCurvesDegraded<unsigned char,unsigned char,unsigned char>(const std::string path,const std::string suffix);
template void LightCurveCollection::writeCurvesDegraded<unsigned char,unsigned short int,unsigned char>(const std::string path,const std::string suffix);
template void LightCurveCollection::writeCurvesDegraded<unsigned short int,unsigned short int,unsigned short int>(const std::string path,const std::string suffix);



LightCurve::LightCurve(int Nsamples){
  this->Nsamples = Nsamples;
  this->t  = (double*) calloc(this->Nsamples,sizeof(double));
  this->m  = (double*) calloc(this->Nsamples,sizeof(double));
  this->dm = (double*) calloc(this->Nsamples,sizeof(double));
}


LightCurve::LightCurve(const LightCurve& other){
  this->Nsamples = other.Nsamples;
  this->t  = (double*) malloc(this->Nsamples*sizeof(double));
  this->m  = (double*) malloc(this->Nsamples*sizeof(double));
  this->dm = (double*) malloc(this->Nsamples*sizeof(double));
  for(int i=0;i<this->Nsamples;i++){
    this->t[i]  = other.t[i];
    this->m[i]  = other.m[i];
    this->dm[i] = other.dm[i];
  }
}

void LightCurve::writeData(const std::string path,const std::string filename){
  std::string full_file = path + filename;
  FILE* fh = fopen(full_file.data(),"w");
  for(int i=0;i<this->Nsamples;i++){
    fprintf(fh,"%11.6e %11.6e %11.6e\n",this->t[i],this->m[i],this->dm[i]);
  }
  fclose(fh);
}

template<typename mType> void LightCurve::writeDegraded(const std::string path,const std::string suffix){
  double m_min,m_max;
  std::string filename = path + "comp_f_" + suffix + ".bin";
  writeQuantity<mType>(filename,this->Nsamples,this->m,m_min,m_max);

  std::string filename_p = path + "comp_p_" + suffix + ".dat";
  double dt = this->t[1] - this->t[0];
  FILE* fh = fopen(filename_p.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"%s\n",typeid(mType).name());
  fprintf(fh,"%f %f\n",m_min,m_max);
  fprintf(fh,"%f\n",dt);
  fclose(fh);
}

template<typename mType,typename tType> void LightCurve::writeDegraded(const std::string path,const std::string suffix){
  double m_min,m_max;
  std::string filename_m = path + "comp_m_" + suffix + ".bin";
  writeQuantity<mType>(filename_m,this->Nsamples,this->m,m_min,m_max);

  double t_min,t_max;
  std::string filename_t = path + "comp_t_" + suffix + ".bin";
  writeQuantity<tType>(filename_t,this->Nsamples,this->t,t_min,t_max);

  std::string filename_p = path + "comp_p_" + suffix + ".dat";
  FILE* fh = fopen(filename_p.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"%s %s\n",typeid(mType).name(),typeid(tType).name());
  fprintf(fh,"%f %f\n",m_min,m_max);
  fprintf(fh,"%f %f\n",t_min,t_max);
  fclose(fh);
}

template<typename mType,typename tType,typename eType> void LightCurve::writeDegraded(const std::string path,const std::string suffix){
  double m_min,m_max;
  std::string filename_m = path + "comp_m_" + suffix + ".bin";
  writeQuantity<mType>(filename_m,this->Nsamples,this->m,m_min,m_max);

  double t_min,t_max;
  std::string filename_t = path + "comp_t_" + suffix + ".bin";
  writeQuantity<tType>(filename_t,this->Nsamples,this->t,t_min,t_max);

  double e_min,e_max;
  std::string filename_e = path + "comp_e_" + suffix + ".bin";
  writeQuantity<eType>(filename_e,this->Nsamples,this->dm,e_min,e_max);
  
  std::string filename_p = path + "comp_p_" + suffix + ".dat";
  FILE* fh = fopen(filename_p.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"'%s' '%s' '%s'\n",typeid(mType).name(),typeid(tType).name(),typeid(eType).name());
  fprintf(fh,"%f %f\n",m_min,m_max);
  fprintf(fh,"%f %f\n",t_min,t_max);
  fprintf(fh,"%f %f\n",e_min,e_max);
  fclose(fh);
}

template<typename qType> void LightCurve::writeQuantity(std::string filename,int Nq,double* q,double& q_min,double& q_max){
  q_min = q[0];
  q_max = q[0];
  for(int i=0;i<Nq;i++){
    if( q[i] > q_max ){
      q_max = q[i];
    }
    if( q[i] < q_min ){
      q_min = q[i];
    }
  }

  std::ofstream out_bin(filename.data(),std::ios::out|std::ios::binary);
  int dum;
  double factor = (q_max-q_min)/((double) std::numeric_limits<qType>::max());
  for(int i=0;i<Nq;i++){
    dum = (int) floor( (q[i]-q_min)/factor );
    out_bin.write((const char*) (&dum),sizeof(qType));
  }
  out_bin.close();
}

template void LightCurve::writeDegraded<unsigned char>(const std::string path,const std::string suffix);
template void LightCurve::writeDegraded<unsigned short int>(const std::string path,const std::string suffix);

template void LightCurve::writeDegraded<unsigned char,unsigned char>(const std::string path,const std::string suffix);
template void LightCurve::writeDegraded<unsigned char,unsigned short int>(const std::string path,const std::string suffix);

template void LightCurve::writeDegraded<unsigned char,unsigned char,unsigned char>(const std::string path,const std::string suffix);
template void LightCurve::writeDegraded<unsigned char,unsigned short int,unsigned char>(const std::string path,const std::string suffix);
template void LightCurve::writeDegraded<unsigned short int,unsigned short int,unsigned short int>(const std::string path,const std::string suffix);

template void LightCurve::writeQuantity<unsigned char>(std::string filename,int Nq,double* q,double& q_min,double& q_max);
template void LightCurve::writeQuantity<unsigned short int>(std::string filename,int Nq,double* q,double& q_min,double& q_max);

