#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <string>
#include <map>
#include "image.hpp"

struct parsSSdisc {
  double mbh;
  double fedd;
  double eta;
};

struct parsParametric {
  double r0;
  double l0;
  double nu;
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
  double width;       // in [10^14 cm]
  double height;      // in [10^14 cm]
  double incl;        // in [degrees]
  double orient;      // in [degrees]
  BaseProfile(double pixSizePhys,double incl,double orient);
  BaseProfile(const BaseProfile& other);
  virtual void generateValues() = 0;
  virtual double getHalfRadius() = 0;
  double sizeParametric(parsParametric pars,double lrest);
  double sizeSS(parsSSdisc pars,double lrest);
protected:
  void makeEven(int& N);
  void createGrid(double* x,double* y);
  void project(int N,double* x,double* y);
  void normalize();
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
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: Gaussian lamp-post ////////////////////////
class GaussianLP : public BaseProfile {
public:
  double sdev; // this is equal to sdev, in [10^14 cm]
  int t0;
  double* f;
  GaussianLP(double pixSizePhys,double sdev,int t0,double* f,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: Gaussian ////////////////////////
class Gaussian : public BaseProfile {
public:
  double sdev; // this is equal to sdev, in [10^14 cm]
  Gaussian(double pixSizePhys,double sdev,double incl,double orient);
  Gaussian(double pixSizePhys,parsParametric pars,double lrest,double incl,double orient);
  Gaussian(double pixSizePhys,parsSSdisc pars,double lrest,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: GaussianHole ////////////////////////
class GaussianHole : public BaseProfile {
public:
  double sdev; // this is equal to sdev, in [10^14 cm]
  double Rin;  // same units as sdev
  GaussianHole(double pixSizePhys,double sdev,double Rin,double incl,double orient);
  GaussianHole(double pixSizePhys,parsParametric pars,double lrest,double Rin,double incl,double orient);
  GaussianHole(double pixSizePhys,parsSSdisc pars,double lrest,double Rin,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: ThermalHole ////////////////////////
class ThermalHole : public BaseProfile {
public:
  double Rin; // in [10^14 cm]
  ThermalHole(double pixSizePhys,double Rin,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: Wavy ////////////////////////
class Wavy : public BaseProfile {
public:
  double a;
  double b;
  int node;
  Wavy(double pixSizePhys,double a,double b,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: Exponential /////////////////
class Exponential : public BaseProfile {
public:
  double sigma;
  Exponential(double pixSizePhys,double sigma,double incl,double orient);
  void generateValues();
  double getHalfRadius();
};

//////////////////////// CLASS DEFINITION: Custom //////////////////////
class Custom : public BaseProfile {
public:
  Custom(double pixSizePhys,const std::string filename,double profPixSizePhys,double incl,double orient);
  void generateValues(){};
  double getHalfRadius(){
    return 0.0;
  };
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

      if( input.shape == "uniform" ){
	return new UniformDisc(input.pixSizePhys,input.pars_parametric,input.lrest,input.incl,input.orient);
      } else if( input.shape == "gaussian" ){
	return new Gaussian(input.pixSizePhys,input.pars_parametric,input.lrest,input.incl,input.orient);
      } else {
	return NULL;
      }

    } else if( input.type == "ss_disc" ){
      
      if( input.shape == "uniform" ){
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

  BaseProfile* createProfileFromMap(std::map<std::string,std::string> input){
    if( input["type"] == "parametric" ){
      parsParametric pars;
      pars.r0 = std::stof(input["r0"]);
      pars.l0 = std::stof(input["l0"]);
      pars.nu = std::stof(input["nu"]);

      if( input["shape"] == "uniform" ){
	return new UniformDisc(std::stof(input["pixSizePhys"]),pars,std::stof(input["lrest"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian" ){
	return new Gaussian(std::stof(input["pixSizePhys"]),pars,std::stof(input["lrest"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else {
	return NULL;
      }

    } else if( input["type"] == "ss_disc" ){
      parsSSdisc pars;
      pars.mbh  = std::stof(input["mbh"]);
      pars.eta  = std::stof(input["eta"]);
      pars.fedd = std::stof(input["fedd"]);
      
      if( input["shape"] == "uniform" ){
	return new UniformDisc(std::stof(input["pixSizePhys"]),pars,std::stof(input["lrest"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian" ){
	return new Gaussian(std::stof(input["pixSizePhys"]),pars,std::stof(input["lrest"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else {
	return NULL;
      }

    } else if( input["type"] == "custom" ){

      return new Custom(std::stof(input["pixSizePhys"]),input["filename"],std::stof(input["profPixSizePhys"]),std::stof(input["incl"]),std::stof(input["orient"]));

    } else {
      return NULL;
    }
  }

private:
  FactoryProfile(){};
};


#endif /* PROFILE_HPP */
