#include "velocities.hpp"

#include <cmath>


velocityComponents::velocityComponents(int N){
  this->N = N;
  this->cmb  = (velocity*) malloc(N*sizeof(velocity));
  this->pec  = (velocity*) malloc(N*sizeof(velocity));
  this->disp = (velocity*) malloc(N*sizeof(velocity));
  this->tot  = (velocity*) malloc(N*sizeof(velocity));
}

void velocityComponents::createVelocitiesK04(int seed,double ra,double dec,double sigma_l,double sigma_s,double sigma_disp,double z_l,double z_s,double D_l,double D_s,double D_ls){
  srand48(seed);

  for(int i=0;i<N;i++){
    velocity vel;

    velCMB(ra,dec,this->cmb[i].v,this->cmb[i].phi,z_l,D_l,D_ls);
    
    this->pec[i].v   = velPec(sigma_l,sigma_s,z_l,z_s,D_l,D_s);
    this->pec[i].phi = getUniform(-180,180);
    
    this->disp[i].v   = velDisp(sigma_disp,z_s,D_s,D_l);
    this->disp[i].phi = getUniform(-180,180);
    
    
    double vtot_ra  = this->pec[i].v*sin(this->pec[i].phi*this->d2r) + this->disp[i].v*sin(this->disp[i].phi*this->d2r) + this->cmb[i].v*sin(this->cmb[i].phi*this->d2r);
    double vtot_dec = this->pec[i].v*cos(this->pec[i].phi*this->d2r) + this->disp[i].v*cos(this->disp[i].phi*this->d2r) + this->cmb[i].v*cos(this->cmb[i].phi*this->d2r);
    this->tot[i].v   = sqrt( pow(vtot_ra,2) + pow(vtot_dec,2) );
    this->tot[i].phi = atan2(vtot_ra,vtot_dec)*this->r2d;
  }
}

void velocityComponents::writeVelocities(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->N;i++){
    fprintf(fh,"%12.4f %7.2f %12.4f %7.2f %12.4f %7.2f %12.4f %7.2f\n",this->tot[i].v,this->tot[i].phi,this->cmb[i].v,this->cmb[i].phi,this->pec[i].v,this->pec[i].phi,this->disp[i].v,this->disp[i].phi);
  }
  fclose(fh);
}


double velocityComponents::getUniform(double a,double b){
  double x = drand48();
  return a+x*(b-a);
}

void velocityComponents::eq2ga(double ra_deg,double dec_deg,double& l_deg,double& b_deg){
  double ra0  = this->ra0_deg*this->d2r;
  double dec0 = this->dec0_deg*this->d2r;
  double l0   = this->l0_deg*this->d2r;
  double ra   = ra_deg*this->d2r;
  double dec  = dec_deg*this->d2r;

  double l = l0 - atan2(cos(dec)*sin(ra-ra0),sin(dec)*cos(dec0)-cos(dec)*sin(dec0)*cos(ra-ra0));
  double b = asin( sin(dec)*sin(dec0)+cos(dec)*cos(dec0)*cos(ra-ra0) );
  l_deg = fmod(l*r2d,360.0);
  b_deg = fmod(b*r2d,360.0);
}

void velocityComponents::ga2eq(double& ra_deg,double& dec_deg,double l_deg,double b_deg){
  double ra0  = this->ra0_deg*this->d2r;
  double dec0 = this->dec0_deg*this->d2r;
  double l0   = this->l0_deg*this->d2r;
  double l    = l_deg*this->d2r;
  double b    = b_deg*this->d2r;

  double ra  = ra0 + atan2(-cos(b)*sin(l-l0),sin(b)*cos(dec0)-cos(b)*sin(dec0)*cos(l-l0));
  double dec = asin( cos(b)*cos(dec0)*cos(l-l0)+sin(b)*sin(dec0) );
  ra_deg  = fmod(ra*this->r2d,360.0);
  dec_deg = fmod(dec*this->r2d,360.0);
}

void velocityComponents::velCMB(double ra_deg,double dec_deg,double& v_cmb,double& phi_cmb,double z_l,double D_l,double D_ls){
  double l_dip     = this->l_dip_deg*this->d2r;
  double b_dip     = this->b_dip_deg*this->d2r;
  

  // convert equatorial to galactic coordinates
  double l_deg;
  double b_deg;
  eq2ga(ra_deg,dec_deg,l_deg,b_deg);
  double l = l_deg*this->d2r;
  double b = (90-b_deg)*this->d2r;
  //std::cout << l << " " << b << std::endl;

  // r is the unit vector in the direction of the observation
  double r_1 = cos(l)*sin(b);
  double r_2 = sin(l)*sin(b);
  double r_3 = cos(b);
  //std::cout << r_1 << " " << r_2 << " " << r_3 << std::endl;

  // vcmb is the vector of the CMB dipole
  double vcmb_1 = this->v_apex*cos(l_dip)*sin(b_dip);
  double vcmb_2 = this->v_apex*sin(l_dip)*sin(b_dip);
  double vcmb_3 = this->v_apex*cos(b_dip);
  //std::cout << vcmb_1 << " " << vcmb_2 << " " << vcmb_3 << std::endl;

  // vnet is the net (transverse) velocity at the direction of observation due to the CMB dipole.
  // It is the vector difference of the velocity of the CMB and its projection at the direction of the observation.
  double a = r_1*vcmb_1 + r_2*vcmb_2 + r_3*vcmb_3;
  double vnet_1 = (vcmb_1 - a*r_1);
  double vnet_2 = (vcmb_2 - a*r_2);
  double vnet_3 = (vcmb_3 - a*r_3);
  //std::cout << vnet_1 << " " << vnet_2 << " " << vnet_3 << std::endl;

  // get vnet back in spherical coordinates
  double norm  = sqrt(vnet_1*vnet_1+vnet_2*vnet_2+vnet_3*vnet_3);
  double b_out = acos(vnet_3/norm)*this->r2d;
  b_out = 90 - b_out;
  double l_out = atan(vnet_2/vnet_1)*this->r2d;
  //std::cout << l_out*d2r << " " << b_out*d2r << std::endl;

  // convert galactic to equatorial coordinates
  double ra_out;
  double dec_out;
  ga2eq(ra_out,dec_out,l_out,b_out);
  //std::cout << norm << "   " << dec_out << "  " << ra_out << std::endl;

  v_cmb = (D_ls/D_l)*norm/(1.0+z_l);
  phi_cmb = atan2(ra_out,dec_out)*this->r2d;
}

double velocityComponents::velPec(double sigma_l,double sigma_s,double z_l,double z_s,double D_l,double D_s){
  double s_term_l = (D_s/D_l)*sigma_l/(1+z_l);
  double s_term_s = sigma_s/(1+z_s);
  double sigma    = sqrt( pow(s_term_l,2) + pow(s_term_s,2) );

  // Applying the Box-Muller transformation
  double u1 = getUniform(0,1);
  double u2 = getUniform(0,1);
  double z1 = sqrt(-2.0 * log(u1)) * cos(this->two_pi * u2);

  return fabs(z1*sigma);
}

double velocityComponents::velDisp(double sigma_disp,double z_l,double D_s,double D_l){
  double e = getUniform(0.8,1.3);
  return sqrt(2)*(D_s/D_l)*e*sigma_disp/(1+z_l);
}



