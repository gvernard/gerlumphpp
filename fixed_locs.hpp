#include <cstdlib>

#include "magnification_map.hpp"

class FixedLocationCollection {
public:
  int Nlocs;
  point* A;
  double* m;
  double* dm;
  int Nx;
  int Ny;
  EffectiveMap* emap;

  FixedLocationCollection(int Nlocs,EffectiveMap* emap);
  FixedLocationCollection(const FixedLocationCollection& other);
  ~FixedLocationCollection(){
    free(A);
    free(m);
    free(dm);
  };


  void setEmap(EffectiveMap* emap);
  void createRandomLocations(int seed);
  void extract();

  void writeLocations(const std::string filename);
  void writeData(const std::string filename);
};
