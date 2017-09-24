#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cstdlib>

#include "image.hpp"
#include "profile.hpp"

class Kernel : public Image{
public:
  int hNx; // half width in pixels of corresponding profile
  int hNy; // half height in pixels of corresponding profile

  Kernel(int map_Nx,int map_Ny,Profile* profile);
  ~Kernel(){
    free(data);
  };
};

#endif /* KERNEL_HPP */
