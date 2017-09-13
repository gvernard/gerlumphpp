#include "map.hpp"
#include "mpd.hpp"

int main(int argc,char* argv[]){

  Map first("12345");


  Mpd* a = first.getFullMpd();


  return 0;
}
