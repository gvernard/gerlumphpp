#include "magnification_map.hpp"
#include "mpd.hpp"
#include "profile.hpp"
#include "kernel.hpp"
#include "light_curve.hpp"
#include "fixed_locs.hpp"

int main(int argc,char* argv[]){
  // This program simply creates to analytic profiles.
  // Can be ran on any machine (no gstar or GERLUMPH data required).


  // Generic options
  double Rein = 500;         // in 10^14 cm
  double pixSizePhys = 1.25; // in 10^14 cm

  // Profile options
  double prof_size_x = 400;  // in 10^14 cm
  double prof_size_y = 400;  // in 10^14 cm

  // A Gaussian profile
  double gauss_sx     = 100;  // in 10^14 cm
  double gauss_sy     = 100;  // in 10^14 cm
  double gauss_incl   = 0;    // in degrees
  double gauss_orient = 0;    // in degrees

  // A Uniform Disc profile
  double disc_radius = 100;  // in 10^14 cm
  double disc_incl   = 0;    // in degrees
  double disc_orient = 0;    // in degrees




  Profile Aprof(1.25,Rein,prof_size_x,prof_size_y);
  Aprof.createGaussian(gauss_sx,gauss_sy,gauss_incl,gauss_orient);
  std::cout << "created Gaussian profile" << std::endl;
  Aprof.writeImageFits("data/gauss.fits",1);
  std::cout << "Gaussian profile written" << std::endl;


  Profile Bprof(1.25,Rein,prof_size_x,prof_size_y);
  Bprof.createUniDisc(disc_radius,disc_incl,disc_orient);
  std::cout << "created Uniform Disc profile" << std::endl;
  Bprof.writeImageFits("data/uniform_disc.fits",1);
  std::cout << "Uniform Disc profile written" << std::endl;




  return 0;
}
