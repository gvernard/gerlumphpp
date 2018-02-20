#include "image.hpp"

void Image::writeImageBIN(const std::string filename,int sampling){
  int nNx = (int) (this->Nx/sampling);
  int nNy = (int) (this->Ny/sampling);

  double dum = 0.0;
  std::ofstream out_bin(filename,std::ios::out|std::ios::binary);
  for(int i=0;i<this->Ny;i+=sampling){
    for(int j=0;j<this->Nx;j+=sampling){
      dum = this->data[i*this->Nx+j];
      out_bin.write((const char*) (&dum),sizeof(double));
    }
  }
  out_bin.close();
}


void Image::writeImageFITS(const std::string filename,int sampling){
  // read, sample, and scale map
  int nNx = (int) (this->Nx/sampling);
  int nNy = (int) (this->Ny/sampling);
  double* output = (double*) calloc(nNx*nNy,sizeof(double));
  long count;

  count = 0;
  for(int i=0;i<this->Ny;i+=sampling){
    for(int j=0;j<this->Nx;j+=sampling){
      output[count] = this->data[i*this->Nx+j];
      count++;
    }
  }

  
  long naxis    = 2;
  long naxes[2] = {(long) nNx,(long) nNy};
  long Ntot = (long) nNx*nNy;
  
  //  std::unique_ptr<CCfits::FITS> pFits(nullptr);
  std::auto_ptr<CCfits::FITS> pFits(0);
  pFits.reset( new CCfits::FITS("!"+filename,FLOAT_IMG,naxis,naxes) );
  
  std::vector<long> extAx(2,(long) nNy);
  CCfits::ExtHDU* imageExt = pFits->addImage("NEW-EXTENSION",FLOAT_IMG,extAx);
  
  //Need Ni and Nj as index counters to flip image
  std::valarray<float> array(Ntot);
  count = 0;
  for(int i=0;i<nNy;i++){
    for(int j=0;j<nNx;j++){
      array[(nNy-1-i)*nNx+j] = (float) (output[count]);
      count++;
    }
  }
  free(output);
  
  long fpixel(1);
  imageExt->write(fpixel,Ntot,array);
  //  pFits->pHDU().addKey("EXPOSURE",13,"Total Exposure Time"); 
  pFits->pHDU().write(fpixel,Ntot,array); 
  //  std::cout << pFits->pHDU() << std::endl;
}

void Image::writeImagePNG(const std::string filename,int sampling){
  // read, sample, and scale map
  int nNx = floor(this->Nx/sampling)-1;
  int nNy = floor(this->Ny/sampling)-1;
  int* colors = (int*) calloc(nNx*nNy,sizeof(int));
  if( this->imageType == "map" ){
    scaleMap(colors,sampling);
  } else if( this->imageType == "profile"){
    scaleProfile(colors,sampling);   
  }
  
  // read rgb values from table file (or select them from a stored list of rgb color tables)
  int* rgb = (int*) calloc(3*256,sizeof(int));
  //  readRGB("rgb.dat",rgb);
  setRGB(rgb);

  // write image
  writeObjectPNG(filename,nNx,nNy,colors,rgb);

  free(colors);
  free(rgb);
}

void Image::scaleMap(int* colors,int sampling){
  double scale_max = 1.6;
  double scale_min = -1.6;
  double scale_fac = 255/(fabs(scale_min) + scale_max);
  double dum,dum2;

  int nNx = floor(this->Nx/sampling)-1;
  int nNy = floor(this->Ny/sampling)-1;

  for(int i=0;i<nNy;i++){
    for(int j=0;j<nNx;j++){
      dum = log10(this->data[i*sampling*this->Nx+j*sampling]);
      if( dum < scale_min ){
	dum = scale_min;
      }
      if( dum > scale_max ){
	dum = scale_max;
      }
      
      dum2 = (dum + fabs(scale_min))*scale_fac;
      colors[i*nNx+j] = (int) round(dum2);
    }
  }
}

void Image::scaleProfile(int* colors,int sampling){
  // Profile is already normalized between 0 and 1
  double scale_max = 0.0;
  double scale_min = 1.0;
  double dum;

  int nNx = floor(this->Nx/sampling)-1;
  int nNy = floor(this->Ny/sampling)-1;

  for(int i=0;i<nNy;i++){
    for(int j=0;j<nNx;j++){
      dum = this->data[i*sampling*this->Nx+j*sampling];
      if( dum < scale_min ){
	scale_min = dum;
      }
      if( dum > scale_max ){
	scale_max = dum;
      }
    }
  }
  double scale_fac = 255/(fabs(scale_min) + scale_max);

  for(int i=0;i<nNy;i++){
    for(int j=0;j<nNx;j++){
      dum = this->data[i*sampling*this->Nx+j*sampling];
      colors[i*nNx+j] = (int) round( (dum + fabs(scale_min))*scale_fac );
    }
  }
}

void Image::readRGB(const std::string filename,int* rgb){
  int r,g,b;
  std::ifstream istr(filename);
  for(int i=0;i<256;i++){
    istr >> r >> g >> b;
    rgb[i*3]     = r;
    rgb[i*3 + 1] = g;
    rgb[i*3 + 2] = b;
  }
}

void Image::setRGB(int* rgb){
  for(int i=0;i<256;i++){
    rgb[i*3]     = this->r[i];
    rgb[i*3 + 1] = this->g[i];
    rgb[i*3 + 2] = this->b[i];
  }
}

void Image::writeObjectPNG(const std::string filename,int width,int height,int* colors,int* rgb){
  FILE* fp            = fopen(filename.data(),"wb");
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  png_infop info_ptr  = png_create_info_struct(png_ptr);
  
  png_init_io(png_ptr,fp);
  png_set_IHDR(png_ptr,info_ptr,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE,PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr,info_ptr);
  
  png_bytep row = (png_bytep) malloc(3 * width * sizeof(png_byte));
  int cindex;
  for(int j=0;j<height;j++) {
    for(int i=0;i<width;i++) {
      cindex = colors[j*width+i];
      
      row[i*3]   = rgb[cindex*3];
      row[i*3+1] = rgb[cindex*3 + 1];
      row[i*3+2] = rgb[cindex*3 + 2];
    }
    png_write_row(png_ptr, row);
  }
  png_write_end(png_ptr, NULL);

  
  fclose(fp);
  png_free_data(png_ptr,info_ptr,PNG_FREE_ALL,-1);
  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
  free(row);
}

