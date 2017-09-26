#ifndef MAGNIFICATION_MAP_HPP
#define MAGNIFICATION_MAP_HPP

#include <cstdlib>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <iterator>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/random.h>
#include <thrust/inner_product.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <cufft.h>

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
  double pixSizePhys; // in units of [10^14 cm]
  bool convolved;

  MagnificationMap(){};
  MagnificationMap(std::string id,double Rein);
  MagnificationMap(const MagnificationMap& other);

  Mpd getFullMpd();
  Mpd getBinnedMpd(int Nbins);
  void convolve(Kernel* kernel,EffectiveMap* emap);
  void writeMapPNG(const std::string filename,int sampling);


private:
  const std::string path = "/lustre/projects/p001_swin/gvernardos/DATABASES/gerlumph_db/";
  //  const std::string path = "/data/users/gvernard/";

  inline void _cudaSafeCall(cudaError err,const char *filename,const int line);
  inline void _cufftSafeCall(cufftResult err,const char* filename,const int line);
  int myfft2d_r2c(int Nx, int Ny, cufftDoubleReal* data, cufftDoubleComplex* Fdata);
  int myfft2d_c2r(int Nx, int Ny, cufftDoubleComplex* Fdata, cufftDoubleReal* data);

  void scaleMap(int* colors,int sampling);
  void readRGB(const std::string filename,int* rgb);
  void writeImage(const std::string filename, int width,int height,int* cols,int* rgb);
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
  Kernel(int map_Nx,int map_Ny,Profile* profile);

  void setKernel(Profile* profile);
};


#endif /* MAGNIFICATION_MAP_HPP */
