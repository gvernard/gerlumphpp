#include "mpd.hpp"

void Mpd::writeMpd(std::string file){
  std::ofstream output(file,std::ios::out);
  for(int i=0;i<this->Nbins;i++){
    output << std::setw(10) << std::left << this->bins[i] << std::setw(10) << std::left << this->counts[i] << std::endl;
  }
  output.close();
}
