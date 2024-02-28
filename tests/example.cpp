#include <string>
#include <vector>

#include "gerlumph.hpp"

using namespace gerlumph;

int main(int argc,char* argv[]){

  // Generic options
  double Rein = 500;        // in 10^14 cm
  std::string map_id = "30090";


  // A Gaussian profile
  double gauss_size   = 100;  // in 10^14 cm
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
  double tmax  = 2000;   // days
  double dt    = 100;    // days
  double v     = 10000;  // km/s
  std::vector<double> t_obs = {10,200,500,800,900,1200,1250,1700};

  // Fixed locations output options
  int Nfixed = 100;
  int fixed_seed = 123;





  std::cout << "Read in a map" << std::endl;
  MagnificationMap amap(map_id,Rein);
  //MagnificationMap bmap("/home/giorgos/myData/","dum_map",Rein);

  std::cout << "Write map images" << std::endl;
  amap.writeImagePNG("sampled_map.png",10);
  amap.writeImageFITS("sampled_map.fits",10);

  
  std::cout << "Get and write map MPDs" << std::endl;
  Mpd a = amap.getFullMpd(); // CPU memory intensive!
  a.writeMpd("full_mpd.dat"); 
  Mpd b = amap.getBinnedMpd(mpd_bins); // CPU memory intensive!
  b.writeMpd("binned_mpd.dat");



  std::cout << "Create a profile and write its image" << std::endl;
  Gaussian aprofile(amap.pixSizePhys,gauss_size,gauss_incl,gauss_orient);
  //UniformDisc aprofoile(amap.pixSizePhys,disc_size,disc_incl,dic_orient);
  aprofile.writeImagePNG("gaussian.png",1);
  aprofile.writeImageFITS("gaussian.fits",1);




  std::cout << "Concolve map with profile" << std::endl;
  Kernel akernel(amap.Nx,amap.Ny,&aprofile);
  EffectiveMap emap(emap_offset,&amap);
  amap.convolve(&akernel,&emap);
  emap.writeImagePNG("convolved_map.png",10);
  emap.writeImageFITS("convolved_map.fits",10);



  
  std::cout << "Get the MPD of the convolved map, and write the map's image" << std::endl;
  emap.writeImagePNG("convolved_map.png",10);
  Mpd c = emap.getBinnedMpd(mpd_bins);
  c.writeMpd("conv_binned_mpd.dat");




  std::cout << "Create a light curve collection with random locations and write them" << std::endl;
  LightCurveCollection acollection(Nlc,&emap);
  acollection.createRandomLocations(lc_seed,lc_length);
  acollection.writeLocations("lc_locs.dat");

  std::cout << "Extract full light curve" << std::endl;
  acollection.extractFull();
  acollection.writeCurves("","lc_full_");

  std::cout << "Copy to a new object (without the extracted data) and extract regularly sampled light curve" << std::endl;
  LightCurveCollection bcollection = acollection;
  bcollection.extractSampled(v,dt,tmax);
  bcollection.writeCurves("","lc_sampled_");

  std::cout << "Copy to a new object (without the extracted data) and extract sampled curve according to observing strategy" << std::endl;
  LightCurveCollection ccollection = acollection;
  ccollection.extractStrategy(v,t_obs);
  ccollection.writeCurves("","lc_strategy_");





  std::cout << "Create a collection of fixed random locations and write output" << std::endl;
  FixedLocationCollection afixed(Nfixed,&emap);
  afixed.createRandomLocations(fixed_seed);
  afixed.extract();
  afixed.writeLocations("fixed_locs_random.dat");
  afixed.writeData("fixed_data_random.dat");

  std::cout << "Create a collection of fixed locations on a grid and write output" << std::endl;
  FixedLocationCollection bfixed(Nfixed,&emap);
  bfixed.createGridLocations();
  bfixed.extract();
  bfixed.writeLocations("fixed_locs_grid.dat");
  bfixed.writeData("fixed_data_grid.dat");





  return 0;
}
