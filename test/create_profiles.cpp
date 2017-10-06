#include <vector>
#include <iostream>
#include "gerlumph.hpp"


int main(int argc,char* argv[]){
  // This program creates and writes two families of uniform discs and Gaussian analytic profiles.
  // Each family can be related by two scaling equations with respect to wavelength: a parametric one and the SS disc one.
  // Can be ran on any machine (no gstar or GERLUMPH data required).


  // Generic options
  double pixSizePhys = 1; // in 10^14 cm
  std::vector<double> lobs{365.49,480.03,622.2,754.06,868.21,992.5}; // observed wavelength in nm
  double z_s = 1.5;

  std::vector<double> lrest;
  for(int i=0;i<lobs.size();i++){
    lrest.push_back( lobs[i]/(1.0+z_s) );
  }

  parsParametric pars{10,102.68,2.33}; // A parametric family of profiles (s0,l0,n)
  //  parsSSdisc pars{0.01,0.3,0.1}; // An SS disc family of profiles (mbh,fedd,eta)




  for(int i=0;i<lrest.size();i++){
    Gaussian gauss(pixSizePhys,pars,lrest[i],0,0);
    //    gauss.writeImagePNG("../data/gauss_"+std::to_string(i)+".png",1);
    gauss.writeImageFITS("../data/gauss_"+std::to_string(i)+".fits",1);

    //    UniformDisc disc(pixSizePhys,pars,lrest[i],0,0);
    //    disc.writeImagePNG("../data/disc_"+std::to_string(i)+".png",1);
  }




  return 0;
}
