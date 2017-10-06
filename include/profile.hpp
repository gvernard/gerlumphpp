#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <string>
#include "image.hpp"

struct parsSSdisc {
  double mbh;
  double fedd;
  double eta;
};

struct parsParametric {
  double s0;
  double l0;
  double n;
};

struct factoryProfilePars {
  double pixSizePhys;
  
  std::string type;
  std::string shape;
  double incl;
  double orient;
  double lrest;

  parsParametric pars_parametric;

  parsSSdisc pars_ssdisc;

  std::string filename;
  double profPixSizePhys;
};

//////////////////////// CLASS DEFINITION: BaseProfile ////////////////////////
class BaseProfile : public Image {
public:
  double pixSizePhys; // in [10^14 cm]
  double incl;
  double orient;
  BaseProfile(double pixSizePhys,double incl,double orient);
  virtual void generateValues() = 0;
protected:
  void makeEven(int& N);
  void normalize();
  double sizeParametric(parsParametric pars,double lrest);
  double sizeSS(parsSSdisc pars,double lrest);
  // check functions here (larger than a map,smaller than a pixel)
  // projection functions here
};

//////////////////////// CLASS DEFINITION: UniformDisc ////////////////////////
class UniformDisc : public BaseProfile {
public:
  double R; // the radius of the disc, in [10^14 cm]
  UniformDisc(double pixSizePhys,double R,double incl,double orient);
  UniformDisc(double pixSizePhys,parsParametric pars,double lrest,double incl,double orient);
  UniformDisc(double pixSizePhys,parsSSdisc pars,double lrest,double incl,double orient);
  void generateValues();
};

//////////////////////// CLASS DEFINITION: Gaussian ////////////////////////
class Gaussian : public BaseProfile {
public:
  double R; // this is equal to sdev, in [10^14 cm]
  Gaussian(double pixSizePhys,double R,double incl,double orient);
  Gaussian(double pixSizePhys,parsParametric pars,double lrest,double incl,double orient);
  Gaussian(double pixSizePhys,parsSSdisc pars,double lrest,double incl,double orient);
  void generateValues();
};

//////////////////////// CLASS DEFINITION: Custom ////////////////////////
class Custom : public BaseProfile {
public:
  Custom(double pixSizePhys,const std::string filename,double profPixSizePhys,double incl,double orient);
  void generateValues(){};
private:
  void interpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
  void binProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
};


//////////////////////// CLASS DEFINITION: FactoryProfile ////////////////////////
class FactoryProfile{ // This is a singleton class.
public:
  FactoryProfile(FactoryProfile const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryProfile const&) = delete;

  static FactoryProfile* getInstance(){
    static FactoryProfile dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseProfile* createProfile(factoryProfilePars input){
    if( input.type == "parametric" ){

      if( input.shape == "uniform_disc" ){
	return new UniformDisc(input.pixSizePhys,input.pars_parametric,input.lrest,input.incl,input.orient);
      } else if( input.shape == "gaussian" ){
	return new Gaussian(input.pixSizePhys,input.pars_parametric,input.lrest,input.incl,input.orient);
      } else {
	return NULL;
      }

    } else if( input.type == "ss_disc" ){
      
      if( input.shape == "uniform_disc" ){
	return new UniformDisc(input.pixSizePhys,input.pars_ssdisc,input.lrest,input.incl,input.orient);
      } else if( input.shape == "gaussian" ){
	return new Gaussian(input.pixSizePhys,input.pars_ssdisc,input.lrest,input.incl,input.orient);
      } else {
	return NULL;
      }

    } else if( input.type == "custom" ){

      return new Custom(input.pixSizePhys,input.filename,input.profPixSizePhys,input.incl,input.orient);

    } else {
      return NULL;
    }
  }

private:
  FactoryProfile(){};
};


#endif /* PROFILE_HPP */
