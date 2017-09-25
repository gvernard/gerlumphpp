#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <string>
#include <math.h>
#include <iostream>

#include <CCfits/CCfits>

#include "image.hpp"

class Profile : public Image {
public:
  double pixSizePhys; // in units of [10^14 cm]
  double width; // in Rein
  double height; // in Rein

  Profile(double pixSizePhys,double Rein,double width,double height);
  Profile(double pixSizePhys,const std::string filename,double profPixSizePhys);
  ~Profile(){
    free(data);
  }

  void createGaussian(double gauss_width,double gauss_height,double incl,double orient);
  void createUniDisc(double radius,double incl,double orient);

private:
  void normalize();
  void interpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
  void binProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);

};

#endif /* PROFILE_HPP */
