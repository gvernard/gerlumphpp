#include "magnification_map.hpp"
#include "mpd.hpp"
#include "profile.hpp"
#include "kernel.hpp"
#include "light_curve.hpp"
#include "fixed_locs.hpp"

int main(int argc,char* argv[]){
  // This program reads a Gaussian "user provided" profile,
  // and interpolates it to have the same resolution as the chosen map (profPixSizePhys > pixSizePhys).
  // It generates an analytical profile with the same characteristics.
  // It convolves both profiles with the same map and compares the resulting light curves.
  // Can be ran ONLY on gstar and requires GERLUMPH data.


  // Generic options
  double Rein = 400;        // in 10^14 cm
  std::string map_id = "18200";


  // Profile options
  double prof_size_x = 400; // in 10^14 cm
  double prof_size_y = 400; // in 10^14 cm

  // A Gaussian profile
  double gauss_sx     = 100;  // in 10^14 cm
  double gauss_sy     = 100;  // in 10^14 cm
  double gauss_incl   = 0;    // in degrees
  double gauss_orient = 0;    // in degrees


  // Effective map output options
  int emap_offset = 1200;   // in pixels


  // Light curve output options
  int Nlc = 10;
  int lc_seed = 123;
  int lc_length = 1600; // in pixels





  // Get magnification map, set emap, set kernel (empty), and light collection (independently of profile)
  MagnificationMap mymap(map_id,Rein);
  EffectiveMap emap(emap_offset,&mymap);
  Kernel mykernel(mymap.Nx,mymap.Ny);
  LightCurveCollection collection(Nlc,&emap);
  collection.createRandomLocations(lc_seed,lc_length);


  // Read custom profile, set kernel, convolve, and extract light curves
  Profile custom_prof(mymap.pixSizePhys,"data/gauss.fits",1.25);
  mykernel.setKernel(&custom_prof);
  mymap.convolve(&mykernel,&emap);
  collection.extractFull();
  collection.writeCurves("custom_full_");


  // Create analytic profile, set kernel, convolve, and extract light curves
  Profile analytic_prof(mymap.pixSizePhys,Rein,prof_size_x,prof_size_y);
  analytic_prof.createGaussian(gauss_sx,gauss_sy,gauss_incl,gauss_orient);
  mykernel.setKernel(&analytic_prof);
  mymap.convolve(&mykernel,&emap);
  collection.extractFull();
  collection.writeCurves("analytic_full_");




  return 0;
}
