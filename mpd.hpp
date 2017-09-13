#ifndef MPD_HPP
#define MPD_HPP

#include <cstdlib>

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

};

#endif /* MPD_HPP */
