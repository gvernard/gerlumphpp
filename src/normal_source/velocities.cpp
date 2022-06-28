#include "velocities.hpp"

#include <cmath>
#include <iostream>

using namespace gerlumph;

velocityComponents::velocityComponents(int N){
  this->N = N;
  this->cmb  = (velocity*) malloc(N*sizeof(velocity));
  this->pec  = (velocity*) malloc(N*sizeof(velocity));
  this->disp = (velocity*) malloc(N*sizeof(velocity));
  this->tot  = (velocity*) malloc(N*sizeof(velocity));

  // Get the dipole location in spherical coordinates.
  // Dipole in equatorial:
  double ra_dip,dec_dip;
  this->ga2eq(ra_dip,dec_dip,this->l_dip_deg,this->b_dip_deg);
  double theta_0 = (90 - dec_dip)*this->d2r;
  double phi_0 = ra_dip*this->d2r;
  // Dipole in galactic:
  //double theta_0 = (90 - this->b_dip_deg)*this->d2r;
  //double phi_0 = this->l_dip_deg*this->d2r;  

  
  // The 3d vector of the CMB velocity
  this->vcmb_x = this->v_apex*cos(phi_0)*sin(theta_0);
  this->vcmb_y = this->v_apex*sin(phi_0)*sin(theta_0);
  this->vcmb_z = this->v_apex*cos(theta_0);
}

void velocityComponents::createVelocitiesK04(int seed,double ra,double dec,double sigma_l,double sigma_s,double sigma_disp,double epsilon,double z_l,double z_s,double D_l,double D_s,double D_ls){
  srand48(seed);

  for(int i=0;i<N;i++){
    velocity vel;

    velCMB(ra,dec,this->cmb[i].v,this->cmb[i].phi,z_l,D_l,D_ls); // the CMB velocity at the RA,DEC of the lens in the RA,DEC reference frame centered on the lens
    //    this->cmb[i].phi -= 90; // reflecting the x-axis of the RA,DEC frame centered on the lens to match the usual right handed cartesian frame
    
    this->pec[i].v   = velPec(sigma_l,sigma_s,z_l,z_s,D_l,D_s);
    this->pec[i].phi = getUniform(-180,180);
    
    this->disp[i].v   = velDisp(sigma_disp,epsilon,z_l,D_s,D_l);
    this->disp[i].phi = getUniform(-180,180);
    
    
    double vtot_x = this->pec[i].v*cos(this->pec[i].phi*this->d2r) + this->disp[i].v*cos(this->disp[i].phi*this->d2r) + this->cmb[i].v*cos(this->cmb[i].phi*this->d2r);
    double vtot_y = this->pec[i].v*sin(this->pec[i].phi*this->d2r) + this->disp[i].v*sin(this->disp[i].phi*this->d2r) + this->cmb[i].v*sin(this->cmb[i].phi*this->d2r);
    this->tot[i].v   = sqrt( pow(vtot_x,2) + pow(vtot_y,2) );
    this->tot[i].phi = atan2(vtot_y,vtot_x)*this->r2d;
  }
}

void velocityComponents::writeVelocities(const std::string filename){
  FILE* fh = fopen(filename.data(),"w");
  for(int i=0;i<this->N;i++){
    this->tot[i].phi  -= 90; // converting to east-of-north
    this->cmb[i].phi  -= 90;
    this->disp[i].phi -= 90;
    this->pec[i].phi  -= 90;
    fprintf(fh,"%12.4f %7.2f %12.4f %7.2f %12.4f %7.2f %12.4f %7.2f\n",this->tot[i].v,this->tot[i].phi,this->cmb[i].v,this->cmb[i].phi,this->disp[i].v,this->disp[i].phi,this->pec[i].v,this->pec[i].phi);
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
  // Convert input ra,dec to spherical coordinates (where ra is aligned with x and dec with the y axes)
  double theta = (90 - dec_deg)*this->d2r;
  double phi   = ra_deg*this->d2r;

  // ur, utheta, and uphi are the local orthogonal unit vectors in the directions of increasing r, theta, and phi
  double ur_x = cos(phi)*sin(theta);
  double ur_y = sin(phi)*sin(theta);
  double ur_z = cos(theta);
  double utheta_x = cos(theta)*cos(phi);
  double utheta_y = cos(theta)*sin(phi);
  double utheta_z = -sin(theta);
  double uphi_x = -sin(phi);
  double uphi_y = cos(phi);
  double uphi_z = 0;

  // v_r is equal to vcmb*r, i.e.  it is the cmb vector projected on the radial unit vector
  double v_r_mag = ur_x*this->vcmb_x + ur_y*this->vcmb_y + ur_z*this->vcmb_z;

  // v_trans is the transverse vector, i.e. it is the vector difference of the CMB vector and its projection at the direction of the observation (radial).
  double v_trans_x = (this->vcmb_x - v_r_mag*ur_x);
  double v_trans_y = (this->vcmb_y - v_r_mag*ur_y);
  double v_trans_z = (this->vcmb_z - v_r_mag*ur_z);
  
  // get magnitude of the transverse velocity vector
  double norm = sqrt(v_trans_x*v_trans_x + v_trans_y*v_trans_y + v_trans_z*v_trans_z);
  v_cmb = (D_ls/D_l)*norm/(1.0+z_l);
  //  v_cmb = norm;
  
  // get the v_r, v_theta, and v_phi coordinates of the transverse vector in the local r,theta,phi reference frame
  // because the transverse vector is by construction on the plane locally tangential to the sphere, vr = 0
  double v_r     = v_trans_x*ur_x     + v_trans_y*ur_y     + v_trans_z*ur_z;
  double v_theta = v_trans_x*utheta_x + v_trans_y*utheta_y + v_trans_z*utheta_z;
  double v_phi   = v_trans_x*uphi_x   + v_trans_y*uphi_y   + v_trans_z*uphi_z;

  // the utheta and uphi are the orthogonal vectors on the tangential plane and they correspond to -x' and -y'
  // hence, the angle of the transverse cmb velocity vector on the plane of the sky defined by x' and y' is:
  phi_cmb = atan2(-v_theta,-v_phi)*this->r2d; // I need -ra because the positive ra axis on the sky is to the left, hence -ra is the coordinate in the usual cartesian frame  
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

double velocityComponents::velDisp(double sigma_disp,double e,double z_l,double D_s,double D_l){
  //  e = getUniform(0.8,1.3);
  return sqrt(2)*(D_s/D_l)*e*sigma_disp/(1+z_l);
}



