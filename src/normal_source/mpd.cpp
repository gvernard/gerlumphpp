#include "mpd.hpp"

using namespace gerlumph;

void Mpd::reset(int Nbins){
  free(this->bins);
  free(this->counts);
  this->Nbins  = Nbins;
  this->bins   = (double*) calloc(Nbins,sizeof(double));
  this->counts = (double*) calloc(Nbins,sizeof(double));
}

void Mpd::writeMpd(const std::string filename){
  /*
  std::ofstream output(file,std::ios::out);
  for(int i=0;i<this->Nbins;i++){
    output << std::setw(10) << std::left << this->bins[i] << std::setw(10) << std::left << this->counts[i] << std::endl;
  }
  output.close();
  */

  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->Nbins;i++){
    fprintf(fh,"%11.6e %11.6e\n",this->bins[i],this->counts[i]);
  }
  fclose(fh);
}
