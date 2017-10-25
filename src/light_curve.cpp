#include "light_curve.hpp"

#include <cstdlib>
#include <cmath>
#include <fstream>

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
template<typename mType> void LightCurveCollection::writeCurvesDegraded(const std::string prefix){
  for(int i=0;i<this->Ncurves;i++){
    std::string suffix = prefix + std::to_string(i);
    this->lightCurves[i].writeDegraded<mType>(suffix);
  }
}
template<typename mType,typename tType> void LightCurveCollection::writeCurvesDegraded(const std::string prefix){
  for(int i=0;i<this->Ncurves;i++){
    std::string suffix = prefix + std::to_string(i);
    this->lightCurves[i].writeDegraded<mType,tType>(suffix);
  }
}
template<typename mType,typename tType,typename eType> void LightCurveCollection::writeCurvesDegraded(const std::string prefix){
  for(int i=0;i<this->Ncurves;i++){
    std::string suffix = prefix + std::to_string(i);
    this->lightCurves[i].writeDegraded<mType,tType,eType>(suffix);
  }
}

void LightCurve::writeData(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Nsamples;i++){
    fprintf(fh,"%11.6e %11.6e %11.6e\n",this->t[i],this->m[i],this->dm[i]);
  }
  fclose(fh);
}

template<typename mType> void LightCurve::writeDegraded(const std::string suffix){
  double m_min,m_max;
  std::string filename = "comp_full_" + suffix + ".bin";
  writeQuantity<mType>(filename,this->Nsamples,this->m,m_min,m_max);

  std::string filename_params = "comp_para_" + suffix + ".dat";
  double dt = this->t[1] - this->t[0];
  FILE* fh = fopen(filename_params.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"%s\n",typeid(mType).name());
  fprintf(fh,"%f %f\n",m_min,m_max);
  fprintf(fh,"%f\n",dt);
  fclose(fh);
}

template<typename mType,typename tType> void LightCurve::writeDegraded(const std::string suffix){
  double m_min,m_max;
  std::string filename_m = "comp_m_" + suffix + ".bin";
  writeQuantity<mType>(filename_m,this->Nsamples,this->m,m_min,m_max);

  double t_min,t_max;
  std::string filename_t = "comp_t_" + suffix + ".bin";
  writeQuantity<tType>(filename_t,this->Nsamples,this->t,t_min,t_max);

  std::string filename_params = "comp_para_" + suffix + ".dat";
  FILE* fh = fopen(filename_params.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"%s %s\n",typeid(mType).name(),typeid(tType).name());
  fprintf(fh,"%f %f\n",m_min,m_max);
  fprintf(fh,"%f %f\n",t_min,t_max);
  fclose(fh);
}

template<typename mType,typename tType,typename eType> void LightCurve::writeDegraded(const std::string suffix){
  double m_min,m_max;
  std::string filename_m = "comp_m_" + suffix + ".bin";
  writeQuantity<mType>(filename_m,this->Nsamples,this->m,m_min,m_max);

  double t_min,t_max;
  std::string filename_t = "comp_t_" + suffix + ".bin";
  writeQuantity<tType>(filename_t,this->Nsamples,this->t,t_min,t_max);

  double e_min,e_max;
  std::string filename_e = "comp_e_" + suffix + ".bin";
  writeQuantity<eType>(filename_e,this->Nsamples,this->dm,e_min,e_max);

  std::string filename_params = "comp_para_" + suffix + ".dat";
  FILE* fh = fopen(filename_params.data(),"w");
  fprintf(fh,"%d\n",this->Nsamples);
  fprintf(fh,"%s %s\n",typeid(mType).name(),typeid(tType).name(),typeid(eType).name());
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
    dum = (int) floor( (q[i]-q_min)*factor);
    out_bin.write((const char*) (&dum),sizeof(qType));
  }
  out_bin.close();
}

