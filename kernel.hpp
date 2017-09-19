#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cstdlib>
#include "profile.hpp"

class Kernel {
public:
  int Nx;  // full kernel width in pixels (same as corresponding map)
  int Ny;  // full kernel height in pixels (same as corresponding map)
  int hNx; // half width in pixels of corresponding profile
  int hNy; // half height in pixels of corresponding profile
  double* data;


  Kernel(int map_Nx,int map_Ny,Profile* profile);
  ~Kernel(){
    free(data);
  };
};

#endif /* KERNEL_HPP */
