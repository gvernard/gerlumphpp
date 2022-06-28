#ifndef FIXED_LOCS_HPP
#define FIXED_LOCS_HPP

#include "magnification_map.hpp"

namespace gerlumph {

  class FixedLocationCollection {
  public:
    int Nlocs;
    std::string type;
    point* A;
    double* m;
    double* dm;
    int Nx;
    int Ny;
    EffectiveMap* emap;

    FixedLocationCollection(int Nlocs,EffectiveMap* emap);
    FixedLocationCollection(int Nlocs,int Nx,int Ny);
    FixedLocationCollection(const FixedLocationCollection& other);
    ~FixedLocationCollection(){
      free(A);
      free(m);
      free(dm);
    };


    void setEmap(EffectiveMap* emap);
    void createRandomLocations(int seed);
    void createGridLocations();
    void extract();

    double checkOverlap(int profRadius);

    void writeLocations(const std::string filename);
    void writeData(const std::string filename);
  };

}
  
#endif /* FIXED_LOCS_HPP */
