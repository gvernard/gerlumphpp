#include <cstdlib>
#include <cmath>
#include <string>

#include "magnification_map.hpp"


class LightCurve {
public:
  int Nsamples;
  double* t;
  double* m;
  double* dm;

  // some functions doing statistics with light curves
};


struct point {
  double x;
  double y;
};

class LightCurveCollection {
public:
  int Ncurves;
  point* A; // initial location of the light curves in pixels
  point* B; // final location of the light curves in pixels
  LightCurve* lightCurves;
  double pixSizePhys;
  int Nx; // width of the effective map from which the light curves will be exracted
  int Ny; // height of the effective map from which the light curves will be exracted

  LightCurveCollection(int Ncurves,EffectiveMap* emap);
  LightCurveCollection(const LightCurveCollection& other);
  ~LightCurveCollection(){
    for(int i=0;i<Ncurves;i++){
      free(lightCurves[i].t);
      free(lightCurves[i].m);
      free(lightCurves[i].dm);
    }
    free(lightCurves);
    free(A);
    free(B);
  }

  void createRandomLocations(int seed,int maxLen);
  void writeLocations(const std::string filename);

  void extractFull(EffectiveMap* emap){};
  void extractSampled(EffectiveMap* emap,double v,double dt,double tmax){};
  void extractStrategy(EffectiveMap* emap,double v,double* t){};


  // sample light curves fully
  // sample light curves according to velocity of the source
  // sample light curves according to observing strategy


};
