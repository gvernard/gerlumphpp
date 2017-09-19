#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <string>
#include <math.h>
#include <iostream>

class Profile {
public:
  double pixSizePhys; // in units of [10^14 cm]
  double width; // in Rein
  double height; // in Rein
  int Nx;
  int Ny;
  double* data;

  Profile(double pixSizePhys,double Rein,double width,double height);
  ~Profile(){
    free(data);
  }

  void readProfile(const std::string filename);
  void createGaussian(double gauss_width,double gauss_height,double incl,double orient);


private:
  void normalize();
};

#endif /* PROFILE_HPP */
