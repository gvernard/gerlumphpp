#ifndef VELOCITIES_HPP
#define VELOCITIES_HPP

#include <cstdlib>
#include <string>

struct velocity{
  double v;
  double phi;
};

class velocityComponents {
public:
  velocity* cmb;
  velocity* pec;
  velocity* disp;
  velocity* tot;
  int N;

  velocityComponents(int N);
  ~velocityComponents(){
    free(cmb);
    free(pec);
    free(disp);
    free(tot);
  };

  void createVelocitiesK04(int seed,double ra,double dec,double sigma_l,double sigma_s,double sigma_disp,double epsilon,double zl,double zs,double Dl,double Ds,double Dls);
  void writeVelocities(const std::string filename);

private:
  const double d2r = 0.017453; // degrees to radians
  const double r2d = 1.0/d2r;  // radians to degrees
  const double ra0_deg  = 192.8595; // ra of the galactic center in deg
  const double dec0_deg = 27.1284;  // dec of the galactic center in deg
  const double l0_deg   = 122.9320; // in deg
  const double two_pi = 2.0*3.14159265358979323846;

  const double v_apex = 387; // km/s
  const double l_dip_deg = 264.4; // in deg
  const double b_dip_deg = 48.4;  // in deg


  void eq2ga(double ra_deg,double dec_deg,double& l_deg,double& b_deg);
  void ga2eq(double& ra_deg,double& dec_deg,double l_deg,double b_deg);
  double getUniform(double a,double b);

  void velCMB(double ra_deg,double dec_deg,double& v_cmb,double& phi_cmb,double z_l,double D_l,double D_ls);
  double velPec(double sigma_l,double sigma_s,double z_l,double z_s,double D_l,double D_s);
  double velDisp(double sigma_disp,double e,double z_l,double D_s,double D_l);
};

#endif /* VELOCITIES_HPP */
