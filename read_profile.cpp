#include "magnification_map.hpp"
#include "mpd.hpp"
#include "profile.hpp"
#include "kernel.hpp"
#include "light_curve.hpp"
#include "fixed_locs.hpp"

int main(int argc,char* argv[]){
  // Generic options
  double Rein = 500;        // in 10^14 cm
  std::string map_id = "18200";

  // Profile options
  double prof_size_x = 400; // in 10^14 cm
  double prof_size_y = 400; // in 10^14 cm

  // A Gaussian profile
  double gauss_sx     = 100;  // in 10^14 cm
  double gauss_sy     = 100;  // in 10^14 cm
  double gauss_incl   = 0;    // in degrees
  double gauss_orient = 0;    // in degrees

  // A Uniform Disc profile
  double disc_radius = 100;  // in 10^14 cm
  double disc_incl   = 0;    // in degrees
  double disc_orient = 0;    // in degrees



  Profile Aprof(1.25,"data/gauss.fits",4);
  std::cout << Aprof.Nx << " " << Aprof.Ny << std::endl;
  Aprof.writeImageFits("data/interpolated.fits",1);


  Profile Bprof(1.25,"data/gauss.fits",0.7);
  std::cout << Bprof.Nx << " " << Bprof.Ny << std::endl;
  Bprof.writeImageFits("data/binned.fits",1);


  return 0;
}
