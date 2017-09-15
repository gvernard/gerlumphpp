#include "map.hpp"
#include "mpd.hpp"
#include "profile.hpp"

int main(int argc,char* argv[]){
  double Rein = 400;       // in 10^14 cm
  double prof_size_x = 20; // in 10^14 cm
  double prof_size_y = 20; // in 10^14 cm


  Map first("67777",Rein);
  //  first.writeMapPNG("testmap.png",1);


  Mpd* a = first.getFullMpd();
  a->writeMpd("full_mpd.dat");
  Mpd* b = first.getBinnedMpd(200);
  b->writeMpd("binned_mpd.dat");


  Profile prof(first.pixSizePhys,Rein,prof_size_x,prof_size_y);
  prof.createGaussian(5,5,0,0);
  double* kernel = prof.setKernel(first.Nx,first.Ny);


  //  first.convolve(kernel);


  return 0;
}
