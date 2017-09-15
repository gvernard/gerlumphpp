#include "mpd.hpp"

void Mpd::writeMpd(const std::string file){
  /*
  std::ofstream output(file,std::ios::out);
  for(int i=0;i<this->Nbins;i++){
    output << std::setw(10) << std::left << this->bins[i] << std::setw(10) << std::left << this->counts[i] << std::endl;
  }
  output.close();
  */

  FILE* fh = fopen(file.data(),"w");
  for(int i=0;i<this->Nbins;i++){
    fprintf(fh,"%11.6e %11.6e\n",this->bins[i],this->counts[i]);
  }
  fclose(fh);

}
