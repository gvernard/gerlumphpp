#ifndef MAP_HPP
#define MAP_HPP

#include <cstdlib>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <fstream>
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
#include <png.h>

#include "mpd.hpp"
#include "profile.hpp"


class LightCurve;


class Map {
public:
  std::string id;
  int Nx;
  int Ny;
  double k;
  double g;
  double s;
  double width; // in Rein
  double height; // in Rein
  double avgmu;
  double avgN;
  double pixSizePhys; // in units of [10^14 cm]
  bool convolved;
  double* data;

  Map(){};
  Map(std::string id,double Rein);
  Map(const Map& other);
  ~Map(){
    free(data);
  };

  Mpd* getFullMpd();
  Mpd* getBinnedMpd(int Nbins);
  void convolve(double* kernel);
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


#endif /* MAP_HPP */
