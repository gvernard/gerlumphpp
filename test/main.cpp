#include "gerlumph.hpp"

int main(int argc,char* argv[]){

  // Generic options
  double Rein = 500;        // in 10^14 cm
  std::string map_id = "18200";


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





  // Read in a map
  MagnificationMap amap(map_id,Rein);

  // Write map images
  amap.writeImagePNG("sampled_map.png",10);
  amap.writeImageFITS("sampled_map.fits",10);

  // Get and write map MPDs
  Mpd a = amap.getFullMpd();
  a.writeMpd("full_mpd.dat");
  Mpd b = amap.getBinnedMpd(mpd_bins);
  b.writeMpd("binned_mpd.dat");





  // Create a profile and write its image
  Gaussian aprofile(amap.pixSizePhys,gauss_size,gauss_incl,gauss_orient);
  aprofile.writeImagePNG("gaussian.png",1);
  aprofile.writeImageFITS("gaussian.fits",1);





  // Concolve map with profile
  Kernel akernel(amap.Nx,amap.Ny,&aprofile);
  EffectiveMap emap(emap_offset,&amap);
  amap.convolve(&akernel,&emap);
  emap.writeImagePNG("convolved_map.png",10);
  emap.writeImageFITS("convolved_map.fits",10);




  // Get the MPD of the convolved map, and write the map's image
  emap.writeImagePNG("convolved_map.png",10);
  Mpd c = emap.getBinnedMpd(mpd_bins);
  c.writeMpd("conv_binned_mpd.dat");





  // Create a light curve collection with random locations and write them
  LightCurveCollection acollection(Nlc,&emap);
  acollection.createRandomLocations(lc_seed,lc_length);
  acollection.writeLocations("lc_locs.dat");

  // Extract full light curve
  acollection.extractFull();
  acollection.writeCurves("lc_full_");
  acollection.writeCurvesDegraded("lc_full_degraded_","byte");

  // Copy to a new object (without the extracted data) and extract regularly sampled light curve
  LightCurveCollection bcollection = acollection;
  bcollection.extractSampled(v,dt,tmax);
  bcollection.writeCurves("lc_sampled_");
  bcollection.writeCurvesDegraded("lc_sampled_degraded_","int16");

  // Copy to a new object (without the extracted data) and extract sampled curve according to observing strategy
  LightCurveCollection ccollection = acollection;
  ccollection.extractStrategy(v,t_obs);
  ccollection.writeCurves("lc_strategy_");
  ccollection.writeCurvesDegraded("lc_strategy_degraded_","int16byte");





  // Create a collection of fixed random locations and write output
  FixedLocationCollection afixed(Nfixed,&emap);
  afixed.createRandomLocations(fixed_seed);
  afixed.extract();
  afixed.writeLocations("fixed_locs_random.dat");
  afixed.writeData("fixed_data_random.dat");

  // Create a collection of fixed locations on a grid and write output
  FixedLocationCollection bfixed(Nfixed,&emap);
  bfixed.createGridLocations();
  bfixed.extract();
  bfixed.writeLocations("fixed_locs_grid.dat");
  bfixed.writeData("fixed_data_grid.dat");





  return 0;
}
