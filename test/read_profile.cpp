#include "magnification_map.hpp"
#include "mpd.hpp"
#include "profile.hpp"
#include "kernel.hpp"
#include "light_curve.hpp"
#include "fixed_locs.hpp"

int main(int argc,char* argv[]){
  // This program reads a "user provided" profile,
  // and interpolates or bins it to the correct pixel resolution (corresponding to the selected map and varying Rein).
  // Can be ran on any machine (no gstar or GERLUMPH data required).


  // Generic options
  double pixSizePhys;
  double profPixSizePhys = 1.25; // in 10^14 cm, can be arbitrary
  int sampling = 1;



  pixSizePhys = 0.666;  // in 10^14 cm
  Profile Aprof(pixSizePhys,"data/gauss.fits",profPixSizePhys);
  Aprof.writeImageFITS("data/interpolated.fits",sampling);
  

  pixSizePhys = 5.1;   // in 10^14 cm
  Profile Bprof(pixSizePhys,"data/gauss.fits",profPixSizePhys);
  Bprof.writeImageFITS("data/binned.fits",sampling);




  return 0;
}
