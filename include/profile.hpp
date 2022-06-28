#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <string>
#include <map>
#include "image.hpp"

namespace gerlumph {

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
    virtual double getHalfRadius() = 0; // for each child profile
    static double getSize(std::map<std::string,std::string> pars,double lrest);
    static double sizeParametric(double r0,double l0,double nu,double lrest);
    static double sizeSS(double mbh,double fedd,double eta,double lrest);
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
    void generateValues();
    double getHalfRadius();
    static double setByHalfRadius(double rhalf);
  };

  //////////////////////// CLASS DEFINITION: Gaussian ////////////////////////
  class Gaussian : public BaseProfile {
  public:
    double sdev; // this is equal to sdev, in [10^14 cm]
    Gaussian(double pixSizePhys,double sdev,double incl,double orient);
    void generateValues();
    double getHalfRadius();
    static double setByHalfRadius(double rhalf);
  };

  //////////////////////// CLASS DEFINITION: GaussianHole ////////////////////////
  class GaussianHole : public BaseProfile {
  public:
    double sdev; // this is equal to sdev, in [10^14 cm]
    double Rin;  // same units as sdev
    GaussianHole(double pixSizePhys,double sdev,double Rin,double incl,double orient);
    void generateValues();
    double getHalfRadius();
    static double setByHalfRadius(double rhalf,double Rin);
  };

  //////////////////////// CLASS DEFINITION: ThermalHole ////////////////////////
  class ThermalHole : public BaseProfile {
  public:
    double Rin; // in [10^14 cm]
    ThermalHole(double pixSizePhys,double Rin,double incl,double orient);
    void generateValues();
    double getHalfRadius();
    static double setByHalfRadius(double rhalf);
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
    static double setByHalfRadius(double rhalf,int node);
  };

  //////////////////////// CLASS DEFINITION: Exponential /////////////////
  class Exponential : public BaseProfile {
  public:
    double sigma;
    Exponential(double pixSizePhys,double sigma,double incl,double orient);
    void generateValues();
    double getHalfRadius();
    static double setByHalfRadius(double rhalf);
  };

  //////////////////////// CLASS DEFINITION: Custom //////////////////////
  class Custom : public BaseProfile {
  public:
    Custom(double pixSizePhys,const std::string filename,double profPixSizePhys,double incl,double orient);
    void generateValues(){};
    double getHalfRadius(){
      // Need to implement the half-ligh radius here
      return 0.0;
    };
    //private:
    //void newInterpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
    //void interpolateProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
    //void binProfile(int Nxx,int Nyy,double* input,double profPixSizePhys);
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

    BaseProfile* createProfileFromPars(std::map<std::string,std::string> input){
      if( input["shape"] == "uniform" ){
	return new UniformDisc(std::stof(input["pixSizePhys"]),std::stof(input["R"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian" ){
	return new Gaussian(std::stof(input["pixSizePhys"]),std::stof(input["sdev"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian_hole" ){
	return new GaussianHole(std::stof(input["pixSizePhys"]),std::stof(input["sdev"]),std::stof(input["Rin"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "thermal_hole" ){
	return new ThermalHole(std::stof(input["pixSizePhys"]),std::stof(input["Rin"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "wavy" ){
	return new Wavy(std::stof(input["pixSizePhys"]),std::stof(input["a"]),std::stof(input["b"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "exponential" ){
	return new Exponential(std::stof(input["pixSizePhys"]),std::stof(input["sigma"]),std::stof(input["incl"]),std::stof(input["orient"]));      
      } else if( input["shape"] == "custom" ){
	return new Custom(std::stof(input["pixSizePhys"]),input["filename"],std::stof(input["profPixSizePhys"]),std::stof(input["incl"]),std::stof(input["orient"]));      
      } else {
	return NULL;
      }
    }

    BaseProfile* createProfileFromHalfRadius(std::map<std::string,std::string> input){
      if( input["shape"] == "uniform" ){
	double R = UniformDisc::setByHalfRadius(std::stof(input["rhalf"]));
	return new UniformDisc(std::stof(input["pixSizePhys"]),R,std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian" ){
	double sdev = Gaussian::setByHalfRadius(std::stof(input["rhalf"]));
	return new Gaussian(std::stof(input["pixSizePhys"]),sdev,std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "gaussian_hole" ){
	double sdev = GaussianHole::setByHalfRadius(std::stof(input["rhalf"]),std::stof(input["Rin"]));
	return new GaussianHole(std::stof(input["pixSizePhys"]),sdev,std::stof(input["Rin"]),std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "thermal_hole" ){
	double Rin = ThermalHole::setByHalfRadius(std::stof(input["rhalf"]));
	return new ThermalHole(std::stof(input["pixSizePhys"]),Rin,std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "wavy" ){
	double b = Wavy::setByHalfRadius(std::stof(input["rhalf"]),3); // for 3 nodes
	return new Wavy(std::stof(input["pixSizePhys"]),std::stof(input["a"]),b,std::stof(input["incl"]),std::stof(input["orient"]));
      } else if( input["shape"] == "exponential" ){
	double sigma = Exponential::setByHalfRadius(std::stof(input["rhalf"]));
	return new Exponential(std::stof(input["pixSizePhys"]),sigma,std::stof(input["incl"]),std::stof(input["orient"]));      
      } else if( input["shape"] == "custom" ){
	return new Custom(std::stof(input["pixSizePhys"]),input["filename"],std::stof(input["profPixSizePhys"]),std::stof(input["incl"]),std::stof(input["orient"]));      
      } else {
	return NULL;
      }      
    }

  private:
    FactoryProfile(){};
  };

}

#endif /* PROFILE_HPP */
