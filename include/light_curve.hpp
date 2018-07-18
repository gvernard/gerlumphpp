#ifndef LIGHT_CURVE_HPP
#define LIGHT_CURVE_HPP

#include <string>
#include <vector>

#include "magnification_map.hpp"

class LightCurve {
public:
  int Nsamples;
  double* t;
  double* m;
  double* dm;


  LightCurve(){};
  LightCurve(int Nsamples);
  LightCurve(const LightCurve& other);
  ~LightCurve(){
    free(t);
    free(m);
    free(dm);
  }

  void writeData(const std::string path,const std::string filename);
  template<typename mType> void writeDegraded(const std::string path,const std::string suffix);
  template<typename mType,typename tType> void writeDegraded(const std::string path,const std::string suffix);
  template<typename mType,typename tType,typename eType> void writeDegraded(const std::string path,const std::string suffix);
  template<typename qType> void writeQuantity(std::string filename,int Nq,double* q,double& q_min,double& q_max);

 // add more functions doing statistics with light curves, maybe get wavelet spectrum, etc
};


class LightCurveCollection {
public:
  int Ncurves;
  std::string type;
  point* A; // initial location of the light curves in pixels
  point* B; // final location of the light curves in pixels
  LightCurve** lightCurves;
  double pixSizePhys;
  int Nx; // width of the effective map from which the light curves will be exracted
  int Ny; // height of the effective map from which the light curves will be exracted
  MagnificationMap* emap;

  LightCurveCollection(){};
  LightCurveCollection(int Ncurves,MagnificationMap* emap);
  LightCurveCollection(const LightCurveCollection& other);
  ~LightCurveCollection(){
    //    for(int i=0;i<this->Ncurves;i++){
      //      delete lightCurves[i];
    //      free(lightCurves[i].t);
    //      free(lightCurves[i].m);
    //      free(lightCurves[i].dm);
    //    }
    //    free(lightCurves);
    delete[] lightCurves;
    free(A);
    free(B);
  }

  void setEmap(MagnificationMap* emap);
  void createRandomLocations(int seed,int maxLen);
  void createOrientedRandomLocations(int seed,int maxLen,double angle);
  void createVelocityLocations(int seed,double tmax,std::vector<double> v,std::vector<double> phi);

  std::vector<int> checkLengthFull();
  std::vector<int> checkLength(double v,double tmax);
  std::vector<int> checkSampling(double v,double dt);
  std::vector<int> checkSampling(double v,std::vector<double> t);

  void extractFull();
  void extractSampled(double v,double dt,double tmax);
  void extractSampled(std::vector<double> v,double dt,double tmax);
  void extractStrategy(double v,std::vector<double> t);
  void extractStrategy(std::vector<double> v,std::vector<double> t);

  void writeLocations(const std::string filename);
  void writeCurves(const std::string path,const std::string suffix);
  template<typename mType> void writeCurvesDegraded(const std::string path,const std::string suffix);
  template<typename mType,typename tType> void writeCurvesDegraded(const std::string path,const std::string suffix);
  template<typename mType,typename tType,typename eType> void writeCurvesDegraded(const std::string path,const std::string suffix);

private:
  const double vfactor = 8.64*1.e-5; // conversion from [km/s] to [10^14 cm / day]
  void sampleLightCurve(int index,std::vector<double> length,double phi);
  void interpolatePlane(double xk,double yk,double& m,double& dm);

  const int full_sampling_lower_limit = 10;
  const int sampling_lower_limit = 3;  
};

#endif /* LIGHT_CURVE_HPP */
