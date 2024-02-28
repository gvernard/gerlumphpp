#include "magnification_map.hpp"

#include <filesystem>

using namespace gerlumph;

MagnificationMap::MagnificationMap():path(MAP_PATH) {
  // do nothing. Just used to define the path.
}

MagnificationMap::MagnificationMap(std::string file_path,std::string file_name,double Rein):path(file_path),id(file_name) {
  std::string full_path = this->path + this->id;
  if( !std::filesystem::exists(full_path+"/map.bin") ){
    throw std::runtime_error("Map file not found at: "+full_path);
  }
  if( !std::filesystem::exists(full_path+"/mapmeta.dat") ){
    throw std::runtime_error("Mapmeta.dat file not found at: "+full_path);
  }
  this->read_map(Rein);
}

MagnificationMap::MagnificationMap(std::string id,double Rein):path(MAP_PATH),id(id) {
  std::string full_path = this->path + this->id;
  if( !std::filesystem::exists(full_path+"/map.bin") ){
    throw std::runtime_error("Map file not found at: "+full_path);
  }
  if( !std::filesystem::exists(full_path+"/mapmeta.dat") ){
    throw std::runtime_error("Mapmeta.dat file not found at: "+full_path);
  }
  this->read_map(Rein);
}

void MagnificationMap::read_map(double Rein){
  std::string full_path = this->path + this->id;
  this->imageType = "map";
  this->convolved = false;
  std::string file;

  // Read map metadata
  file = full_path + "/mapmeta.dat";
  std::ifstream myfile(file.c_str());
  myfile >> this->avgmu >> this->avgN;
  myfile >> this->Nx;
  myfile >> this->width;
  myfile >> this->k >> this->g >> this->s;
  myfile.close();
  this->Ny = this->Nx;
  this->height = this->width;
  this->pixSizePhys = Rein*this->width/this->Nx; // in units of [10^14 cm]
  this->mu_th = 1.0/(pow(1.0-this->k,2)-pow(this->g,2));

  // Read map data
  file = full_path + "/map.bin";
  FILE* ptr_myfile = fopen(file.data(),"rb");
  int* imap = (int*) calloc(this->Nx*this->Ny,sizeof(int));
  size_t dum = fread(imap,sizeof(int),this->Nx*this->Ny,ptr_myfile);
  fclose(ptr_myfile);
  
  double factor = fabs(this->avgmu/this->avgN);
  double muth   = fabs( 1.0/(pow(1.0-this->k,2)-pow(this->g,2)) );
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = (double) (imap[i]*factor/muth);
  }
  free(imap);
}


MagnificationMap::MagnificationMap(const MagnificationMap& other):path(other.path){
  this->imageType = other.imageType;
  this->id     = other.id;
  this->k      = other.k;
  this->g      = other.g;
  this->s      = other.s;
  this->Nx     = other.Nx;
  this->Ny     = other.Ny;
  this->width  = other.width;
  this->height = other.height;
  this->avgmu  = other.avgmu;
  this->avgN   = other.avgN;
  this->mu_th  = other.mu_th;
  this->pixSizePhys = other.pixSizePhys; // in units of [10^14 cm]
  this->convolved   = other.convolved;

  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  for(long i=0;i<this->Nx*this->Ny;i++){
    this->data[i] = other.data[i];
  }
}


double MagnificationMap::getPixSizePhys(std::string id,double Rein){
  // Read map metadata
  std::string file = MAP_PATH + id + "/mapmeta.dat";
  std::ifstream myfile(file.c_str());
  double avgmu,avgN,Nx,width;
  myfile >> avgmu >> avgN;
  myfile >> Nx; // it is an int, but it doesn't matter here
  myfile >> width;
  return Rein*width/Nx; // in units of [10^14 cm]
}


EffectiveMap::EffectiveMap(int offset,MagnificationMap* map){
  this->top    = offset;
  this->bottom = offset;
  this->left   = offset;
  this->right  = offset;

  this->imageType = map->imageType;
  this->pixSizePhys = map->pixSizePhys;
  this->Nx = map->Nx - 2*offset;
  this->Ny = map->Ny - 2*offset;
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width = map->width*this->Nx/map->Nx;
  this->height = map->height*this->Ny/map->Ny;

  this->k = map->k;
  this->g = map->g;
  this->s = map->s;
  this->avgmu = map->avgmu;
  this->avgN = map->avgN;

  this->convolved = true;
}

EffectiveMap::EffectiveMap(double d_offset,MagnificationMap* map){
  this->pixSizePhys = map->pixSizePhys;
  int offset = ceil(d_offset/this->pixSizePhys);

  this->top    = offset;
  this->bottom = offset;
  this->left   = offset;
  this->right  = offset;

  this->imageType = map->imageType;
  this->Nx = map->Nx - 2*offset;
  this->Ny = map->Ny - 2*offset;
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width = map->width*this->Nx/map->Nx;
  this->height = map->height*this->Ny/map->Ny;

  this->k = map->k;
  this->g = map->g;
  this->s = map->s;
  this->avgmu = map->avgmu;
  this->avgN = map->avgN;

  this->convolved = true;
}

EffectiveMap::EffectiveMap(int top,int bottom,int left,int right,MagnificationMap* map){
  this->top    = top;
  this->bottom = bottom;
  this->left   = left;
  this->right  = right;

  this->imageType = map->imageType;
  this->pixSizePhys = map->pixSizePhys;
  this->Nx = map->Nx - left - right;
  this->Ny = map->Ny - top - bottom;
  this->data = (double*) calloc(this->Nx*this->Ny,sizeof(double));
  this->width = map->width*this->Nx/map->Nx;
  this->height = map->height*this->Ny/map->Ny;

  this->k = map->k;
  this->g = map->g;
  this->s = map->s;
  this->avgmu = map->avgmu;
  this->avgN = map->avgN;

  this->convolved = true;
}





Kernel::Kernel(int map_Nx,int map_Ny){
  this->imageType = "map";
  this->Nx   = map_Nx;
  this->Ny   = map_Ny;
  this->data = (double*) calloc(map_Nx*map_Ny,sizeof(double));
}

Kernel::Kernel(int map_Nx,int map_Ny,BaseProfile* profile){
  this->imageType = "map";
  this->Nx   = map_Nx;
  this->Ny   = map_Ny;
  this->data = (double*) calloc(map_Nx*map_Ny,sizeof(double));
  setKernel(profile);
}

void Kernel::setKernel(BaseProfile* profile){
  this->hNx = profile->Nx/2;
  this->hNy = profile->Ny/2;

  int fNx = profile->Nx;   // full profile width in pixels
  //  int fNy = profile->Ny;   // full profile height in pixels
  int hNx = profile->Nx/2; // half profile width in pixels
  int hNy = profile->Ny/2; // half profile height in pixels
  int map_Nx = this->Nx;
  int map_Ny = this->Ny;

  free(data);
  this->data = (double*) calloc(map_Nx*map_Ny,sizeof(double));

  for(int i=0;i<hNy;i++){
    for(int j=0;j<hNx;j++){
      this->data[i*map_Nx+j]                                = profile->data[hNy*fNx+hNx+i*fNx+j];
      this->data[map_Nx-hNx+i*map_Nx+j]                     = profile->data[hNy*fNx+i*fNx+j];
      this->data[map_Nx*(map_Ny-hNy)+i*map_Nx+j]            = profile->data[hNx+fNx*i+j];
      this->data[map_Nx*(map_Ny-hNy)+map_Nx-hNx+i*map_Nx+j] = profile->data[fNx*i+j];
    }
  }
}
