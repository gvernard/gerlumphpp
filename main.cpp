#include "magnification_map.hpp"
#include "mpd.hpp"
#include "profile.hpp"
#include "kernel.hpp"
#include "light_curve.hpp"

int main(int argc,char* argv[]){
  double Rein = 400;        // in 10^14 cm
  double prof_size_x = 1200; // in 10^14 cm
  double prof_size_y = 1200; // in 10^14 cm
  double gauss_sx  = 400;    // in 10^14 cm
  double gauss_sy  = 400;    // in 10^14 cm
  double gauss_incl = 0;    // in degrees
  double gauss_orient = 0;  // in degrees

  int emap_offset = 1200;   // in pixels
  int lc_length = 400;     // in pixels




  MagnificationMap amap("67777",Rein);

  //  amap.writeMapPNG("full_map.png",1);
  //  std::cout << "full map written" << std::endl;
  //  amap.writeMapPNG("sampled_map.png",10);
  //  std::cout << "map sampled" << std::endl;

  //  Mpd a = amap.getFullMpd();
  //  a.writeMpd("full_mpd.dat");
  //  std::cout << "full mpd written" << std::endl;


  Mpd b = amap.getBinnedMpd(200);
  b.writeMpd("binned_mpd.dat");
  std::cout << "binned mpd written" << std::endl;


  Profile aprofile(amap.pixSizePhys,Rein,prof_size_x,prof_size_y);
  aprofile.createGaussian(gauss_sx,gauss_sy,gauss_incl,gauss_orient);
  std::cout << "created Gaussian profile" << std::endl;


  Kernel akernel(amap.Nx,amap.Ny,&aprofile);
  std::cout << "created kernel based on the profile and current map size" << std::endl;


  EffectiveMap emap(emap_offset,&amap);
  std::cout << "created effective map based on current map size" << std::endl;


  amap.convolve(&akernel,&emap);
  std::cout << "convolved map with kernel and profile" << std::endl;


  //  emap.writeMapPNG("convolved_map.png",10);
  //  std::cout << "convolved map written" << std::endl;


  //  Mpd c = emap.getBinnedMpd(200);
  //  c.writeMpd("conv_binned_mpd.dat");
  //  std::cout << "convolved and binned mpd written" << std::endl;


  LightCurveCollection acollection(10,&emap);
  acollection.createRandomLocations(123,lc_length);
  acollection.writeLocations("lc_locs.dat");

  return 0;
}
