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



  // Effective map output options
  int emap_offset = 1200;   // in pixels

  // MPD output options
  int mpd_bins = 200;

  // Light curve output options
  int Nlc = 10;
  int lc_seed = 123;
  int lc_length = 1600; // in pixels
  double tmax = 2000;   // days
  double dt   = 100;    // days
  double v    = 10000;  // km/s

  // Fixed locations output options
  int Nfixed = 100;
  int fixed_seed = 123;






  MagnificationMap amap(map_id,Rein);

  //  amap.writeImagePNG("full_map.png",1);
  //  std::cout << "full map written" << std::endl;
  //  amap.writeImagePNG("sampled_map.png",10);
  //  std::cout << "map sampled" << std::endl;

  //  Mpd a = amap.getFullMpd();
  //  a.writeMpd("full_mpd.dat");
  //  std::cout << "full mpd written" << std::endl;


  Mpd b = amap.getBinnedMpd(mpd_bins);
  b.writeMpd("binned_mpd.dat");
  std::cout << "binned mpd written" << std::endl;


  Profile aprofile(amap.pixSizePhys,Rein,prof_size_x,prof_size_y);
  aprofile.createGaussian(gauss_sx,gauss_sy,gauss_incl,gauss_orient);
  std::cout << "created Gaussian profile" << std::endl;
  aprofile.writeImagePNG("gaussian.png",1);
  std::cout << "Gaussian profile written" << std::endl;


  Profile bprofile(amap.pixSizePhys,Rein,prof_size_x,prof_size_y);
  bprofile.createUniDisc(disc_radius,disc_incl,disc_orient);
  std::cout << "created Uniform Disc profile" << std::endl;
  bprofile.writeImagePNG("uniform_disc.png",1);
  std::cout << "Uniform Disc profile written" << std::endl;




  Kernel akernel(amap.Nx,amap.Ny,&aprofile);
  std::cout << "created kernel based on the profile and current map size" << std::endl;


  EffectiveMap emap(emap_offset,&amap);
  std::cout << "created effective map based on current map size" << std::endl;


  amap.convolve(&akernel,&emap);
  std::cout << "convolved map with kernel and profile" << std::endl;


  emap.writeImagePNG("convolved_map.png",10);
  std::cout << "convolved map written" << std::endl;


  //  Mpd c = emap.getBinnedMpd(mpd_bins);
  //  c.writeMpd("conv_binned_mpd.dat");
  //  std::cout << "convolved and binned mpd written" << std::endl;




  LightCurveCollection acollection(Nlc,&emap);
  acollection.createRandomLocations(lc_seed,lc_length);
  std::cout << "light curve locations created" << std::endl;
  acollection.writeLocations("lc_locs.dat");
  std::cout << "light curve locations written" << std::endl;
  acollection.extractFull();
  std::cout << "light curve data extracted" << std::endl;
  acollection.writeCurves("lcdata_full_");
  std::cout << "light curve data written" << std::endl;




  FixedLocationCollection fixed(Nfixed,&emap);
  fixed.createRandomLocations(fixed_seed);
  fixed.extract();
  fixed.writeLocations("fixed_locs.dat");
  fixed.writeData("fixed_data.dat");



  /*
  LightCurveCollection bcollection = acollection;
  std::cout << "light curve collection copied" << std::endl;
  bcollection.extractSampled(v,dt,tmax);
  std::cout << "light curves extracted" << std::endl;
  for(int i=0;i<bcollection.Ncurves;i++){
    std::string fname = "lcdata_sampled_" + std::to_string(i) + ".dat";
    bcollection.lightCurves[i].writeData(fname);
  }
  std::cout << "light curves written" << std::endl;
  */

  /*
  double Lmax = 8.64*1.e-5*v*tmax/emap.pixSizePhys;
  LightCurveCollection bcollection(Nlc,&emap);
  bcollection.createRandomLocations(lc_seed,floor(Lmax));
  std::cout << "light curve positions created" << std::endl;
  bcollection.extractSampled(v,dt,tmax);
  std::cout << "light curves extracted" << std::endl;
  for(int i=0;i<bcollection.Ncurves;i++){
    std::string fname = "lcdata_sampled_" + std::to_string(i) + ".dat";
    bcollection.lightCurves[i].writeData(fname);
  }
  std::cout << "light curves written" << std::endl;
  */

  /*
  std::vector<double> t = {100,120,140,160,180,200,500,1000,1500,1800};
  double Lmax2 = 8.64*1.e-5*v*t[9]/emap.pixSizePhys;
  LightCurveCollection ccollection(Nlc,&emap);
  ccollection.createRandomLocations(lc_seed,floor(Lmax2));
  std::cout << "light curve positions created" << std::endl;
  ccollection.extractStrategy(v,t);
  std::cout << "light curves extracted" << std::endl;
  for(int i=0;i<ccollection.Ncurves;i++){
    std::string fname = "lcdata_strategy_" + std::to_string(i) + ".dat";
    ccollection.lightCurves[i].writeData(fname);
  }
  std::cout << "light curves written" << std::endl;
  */






  return 0;
}
