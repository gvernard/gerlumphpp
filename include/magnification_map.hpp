#ifndef MAGNIFICATION_MAP_HPP
#define MAGNIFICATION_MAP_HPP

#ifdef MAP_PATH
  #define PATH MAP_PATH
#else
//  #define PATH "/nfs/cluster/cosmic/gvernardos/DATABASES/gerlumph_db/"
  #define PATH "/lustre/projects/p001_swin/gvernardos/DATABASES/gerlumph_db/"
#endif


#include <cstdlib>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <iterator>

#include "image.hpp"
#include "mpd.hpp"
#include "profile.hpp"

class EffectiveMap;
class Kernel;

struct point {
  // required in light_curve.hpp and fixed_locs.hpp where this file is included as a header
  double x;
  double y;
};

class MagnificationMap : public Image {
public:
  std::string id;
  double k;
  double g;
  double s;
  double width; // in Rein
  double height; // in Rein
  double avgmu;
  double avgN;
  double mu_th;
  double pixSizePhys; // in units of [10^14 cm]
  bool convolved;

  MagnificationMap(){};
  MagnificationMap(std::string id,double Rein);
  MagnificationMap(const MagnificationMap& other);

  Mpd getFullMpd();
  Mpd getBinnedMpd(int Nbins);
  void convolve(Kernel* kernel,EffectiveMap* emap);


private:
  const std::string path = PATH;
  //  const std::string path = "/data/users/gvernard/";
};


class EffectiveMap : public MagnificationMap {
public:
  int top;    // top margin in pixels
  int bottom; // bottom margin in pixels
  int left;   // left margin in pixels
  int right;  // right margin in pixels


  EffectiveMap(int offset,MagnificationMap* map);
  EffectiveMap(double d_offset,MagnificationMap* map);
  EffectiveMap(int top,int bottom,int left,int right,MagnificationMap* map);
};


class Kernel : public Image{
public:
  int hNx; // half width in pixels of corresponding profile
  int hNy; // half height in pixels of corresponding profile

  Kernel(int map_Nx,int map_Ny);
  Kernel(int map_Nx,int map_Ny,BaseProfile* profile);

  void setKernel(BaseProfile* profile);
};


#endif /* MAGNIFICATION_MAP_HPP */
