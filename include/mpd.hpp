#ifndef MPD_HPP
#define MPD_HPP

#include <cstdlib>
#include <string>
#include <fstream>
#include <iomanip>

namespace gerlumph {

  class Mpd {
  public:
    int Nbins;
    double* bins;
    double* counts;

    Mpd(int Nbins){
      this->Nbins  = Nbins;
      this->bins   = (double*) calloc(Nbins,sizeof(double));
      this->counts = (double*) calloc(Nbins,sizeof(double));
    }
    ~Mpd(){
      free(bins);
      free(counts);
    }

    void reset(int Nbins);
    void writeMpd(const std::string filename);
  };

}
  
#endif /* MPD_HPP */
