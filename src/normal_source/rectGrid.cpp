#include "rectGrid.hpp"
#include "fitsInterface.hpp"

#include <cmath>
#include <algorithm> // for transform

#include <iostream>


// Definition of static variables:
// 1st derivatives
const std::vector<int>    Constants::derivative_1_forward_1_index({0,1});
const std::vector<double> Constants::derivative_1_forward_1_coeff({-1,1});
const std::vector<int>    Constants::derivative_1_forward_2_index({0,1,2});
const std::vector<double> Constants::derivative_1_forward_2_coeff({-1.5,2.0,-0.5});
const std::vector<int>    Constants::derivative_1_backward_1_index({-1,0});
const std::vector<double> Constants::derivative_1_backward_1_coeff({-1,1});
const std::vector<int>    Constants::derivative_1_backward_2_index({-2,-1,0});
const std::vector<double> Constants::derivative_1_backward_2_coeff({0.5,-2,1.5});
const std::vector<int>    Constants::derivative_1_central_2_index({-1,0,1});
const std::vector<double> Constants::derivative_1_central_2_coeff({-0.5,0.0,0.5});
// 2nd derivatives
const std::vector<int>    Constants::derivative_2_forward_1_index({0,1,2});
const std::vector<double> Constants::derivative_2_forward_1_coeff({1,-2,1});  
const std::vector<int>    Constants::derivative_2_forward_2_index({0,1,2,3});
const std::vector<double> Constants::derivative_2_forward_2_coeff({2,-5,4,-1});
const std::vector<int>    Constants::derivative_2_backward_1_index({-2,-1,0});
const std::vector<double> Constants::derivative_2_backward_1_coeff({1,-2,1});  
const std::vector<int>    Constants::derivative_2_backward_2_index({-3,-2,-1,0});
const std::vector<double> Constants::derivative_2_backward_2_coeff({-1,4,-5,2});
const std::vector<int>    Constants::derivative_2_central_2_index({-1,0,1});
const std::vector<double> Constants::derivative_2_central_2_coeff({1,-2,1});


RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::map<std::string,std::string> options){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax,options);
}

RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,const std::string filepath,std::map<std::string,std::string> options){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax,options);
  FitsInterface::readFits(this->Ny,this->Nx,this->z,filepath);
}

void RectGrid::common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::map<std::string,std::string> options){
  this->Nx = Nx;
  this->Ny = Ny;
  this->width  = xmax - xmin;
  this->height = ymax - ymin;
  this->xmin = xmin;
  this->xmax = xmax;
  this->ymin = ymin;
  this->ymax = ymax;

  this->step_x = this->width/this->Nx;
  this->step_y = this->height/this->Ny;

  for(std::map<std::string,std::string>::iterator it=options.begin();it!=options.end();it++){
    this->options[it->first] = it->second;
  }
  this->set_interp(options["interp"]);
  
  this->Nz = this->Nx*this->Ny;
  this->z  = (double*) calloc(this->Nx*this->Ny,sizeof(double));

  this->center_y = (double*) calloc(this->Ny,sizeof(double));
  this->bound_y  = (double*) calloc(this->Ny+1,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    this->center_y[i] = ymin + (this->step_y/2.0) + i*this->step_y;
    this->bound_y[i]  = ymin + i*this->step_y;
  }  
  this->bound_y[this->Ny] = ymax;
  
  this->center_x = (double*) calloc(this->Nx,sizeof(double));
  this->bound_x  = (double*) calloc(this->Nx+1,sizeof(double));
  for(int j=0;j<this->Nx;j++){
    this->center_x[j] = xmin + (this->step_x/2.0) + j*this->step_x;
    this->bound_x[j]  = xmin + j*this->step_x;
  }
  this->bound_x[this->Nx] = xmax;
}

void RectGrid::set_interp(std::string interp){
  if( interp == "bilinear" ){
    this->interp2d = &RectGrid::interp2d_bilinear;
  } else if( interp == "bicubic" ){
    this->interp2d = &RectGrid::interp2d_bicubic;
  } else {
    this->interp2d = &RectGrid::interp2d_nearest;
  }
}


RectGrid::RectGrid(const RectGrid& grid){
  this->common_constructor(grid.Nx,grid.Ny,grid.xmin,grid.xmax,grid.ymin,grid.ymax,grid.options);
  for(int k=0;k<grid.Nz;k++){
    this->z[k] = grid.z[k];
  }
  if( grid.zx != NULL ){
    for(int k=0;k<grid.Nz;k++){
      this->zx = (double*) calloc(this->Nz,sizeof(double));
      this->zx[k] = grid.zx[k];
    }
  }
  if( grid.zy != NULL ){
    for(int k=0;k<grid.Nz;k++){
      this->zy = (double*) calloc(this->Nz,sizeof(double));
      this->zy[k] = grid.zy[k];
    }
  }
  if( grid.zxy != NULL ){
    for(int k=0;k<grid.Nz;k++){
      this->zxy = (double*) calloc(this->Nz,sizeof(double));
      this->zxy[k] = grid.zxy[k];
    }
  }
  if( grid.zxx != NULL ){
    for(int k=0;k<grid.Nz;k++){
      this->zxx = (double*) calloc(this->Nz,sizeof(double));
      this->zxx[k] = grid.zxx[k];
    }
  }
  if( grid.zyy != NULL ){
    for(int k=0;k<grid.Nz;k++){
      this->zyy = (double*) calloc(this->Nz,sizeof(double));
      this->zyy[k] = grid.zyy[k];
    }
  }  
}

RectGrid::~RectGrid(){
  free(center_x);
  free(center_y);
  free(z);
  free(zx);
  free(zy);
  free(zxy);
  free(zxx);
  free(zyy);
  free(bound_x);
  free(bound_y);
}

RectGrid RectGrid::embeddedNewGrid(int new_Nx,int new_Ny,std::string mode){
  double new_xmin = this->bound_x[0];
  double new_xmax = this->bound_x[this->Nx];
  double new_ymin = this->bound_y[0];
  double new_ymax = this->bound_y[this->Ny];
  RectGrid new_grid(new_Nx,new_Ny,new_xmin,new_xmax,new_ymin,new_ymax,this->options);

  if( this->z == NULL ){
    // std::cout << "z is NULL" << std::endl;
    return new_grid;
  } else if( new_Nx > this->Nx && new_Ny > this->Ny ){
    new_grid.z = (double*) calloc(new_grid.Nz,sizeof(double));
    for(int i=0;i<new_Ny;i++){
      for(int j=0;j<new_Nx;j++){
	double x = new_grid.center_x[j];
	double y = new_grid.center_y[i];
	if( this->point_between_pixel_centers(x,y,0) ){
	  new_grid.z[i*new_grid.Nx+j] = (this->*interp2d)(x,y,this->z);
	} else {
	  int i0,j0;
	  match_point_to_pixel(x,y,i0,j0);
	  new_grid.z[i*new_grid.Nx+j] = this->z[i0*this->Nx+j0];
	}
      }
    }
    return new_grid;
  } else if( new_Nx < this->Nx && new_Ny < this->Ny ){
    new_grid.z = (double*) calloc(new_grid.Nz,sizeof(double));
    if( mode == "interp" ){
      for(int i=0;i<new_Ny;i++){
	for(int j=0;j<new_Nx;j++){
	  double x = new_grid.center_x[j];
	  double y = new_grid.center_y[i];
	  new_grid.z[i*new_grid.Nx+j] = (this->*interp2d)(x,y,this->z);
	}
      }
      return new_grid;
    } else if( mode == "additive" ){
      int i0,j0;
      for(int i=0;i<this->Ny;i++){
	for(int j=0;j<this->Nx;j++){
	  new_grid.match_point_to_pixel(this->center_x[j],this->center_y[i],i0,j0);
	  new_grid.z[i0*new_grid.Nx + j0] += this->z[i*this->Nx + j];
	}
      }
      return new_grid;
    } else if( mode == "integrate" ){
      int* counts = (int*) calloc(new_grid.Nz,sizeof(int));
      int i0,j0;
      for(int i=0;i<this->Ny;i++){
	for(int j=0;j<this->Nx;j++){
	  new_grid.match_point_to_pixel(this->center_x[j],this->center_y[i],i0,j0);
	  new_grid.z[i0*new_grid.Nx + j0] += this->z[i*this->Nx + j];
	  counts[i0*new_grid.Nx + j0] += 1;
	}
      }
      double new_pix_area = new_grid.step_x*new_grid.step_y;
      for(int i=0;i<new_grid.Nz;i++){
	new_grid.z[i] = new_pix_area*new_grid.z[i]/counts[i];
      }
      free(counts);
      return new_grid;
    } else {
      //      std::cout << "something wrong" << std::endl;
      //throw exception
      return new_grid;
    }

  } else {
    // std::cout << "dimensions wrong" << std::endl;
    // throw exception
    return new_grid;	  
  }
}

bool RectGrid::point_in_grid(double x,double y){
  if( x < this->xmin || this->xmax < x || y < this->ymin || this->ymax < y ){
    return false;
  } else {
    return true;
  }
}

bool RectGrid::point_between_pixel_centers(double x,double y,int boundary_size){
  if( x < this->center_x[0] || this->center_x[this->Nx-1-boundary_size] < x || y < this->center_y[0] || this->center_y[this->Ny-1-boundary_size] < y ){
    return false;
  } else {
    return true;
  }
}

void RectGrid::match_point_to_pixel(double x,double y,int& i0,int& j0){
  j0 = (int) floor( (x-this->xmin)/this->step_x );
  i0 = (int) floor( (y-this->ymin)/this->step_y );
}

bool RectGrid::match_point_to_closest_4(double x,double y,int* i,int* j){
  if( !this->point_between_pixel_centers(x,y,0) ){
    return false;
  }
  int i0,j0;
  match_point_to_pixel(x,y,i0,j0);
  double dy = y - this->center_y[i0];
  double dx = x - this->center_x[j0];
  int f1,f2;
  if( dx<0 ){
    f1 = -1;
    f2 = 0;
  } else {
    f1 = 0;
    f2 = 1;
  }
  j[0] = j0 + f1;
  j[1] = j0 + f2;
  j[2] = j0 + f1;
  j[3] = j0 + f2;
  if( dy<0 ){
    f1 = 0;
    f2 = -1;
  } else {
    f1 = 1;
    f2 = 0;
  }
  i[0] = i0 + f1;
  i[1] = i0 + f1;
  i[2] = i0 + f2;
  i[3] = i0 + f2;  
  return true;
}

bool RectGrid::match_point_to_closest_16(double x,double y,int* i,int* j){
  if( !this->point_between_pixel_centers(x,y,1) ){
    return false;
  }
  int i0,j0;
  match_point_to_pixel(x,y,i0,j0);
  double dy = y - this->center_y[i0];
  double dx = x - this->center_x[j0];
  int f1,f2,f3,f4;
  if( dx<0 ){
    f1 = -2;
    f2 = -1;
    f3 = 0;
    f4 = 1;
  } else {
    f1 = -1;
    f2 = 0;
    f3 = 1;
    f4 = 2;
  }
  j[0]  = j0 + f1;
  j[1]  = j0 + f2;
  j[2]  = j0 + f3;
  j[3]  = j0 + f4;
  j[4]  = j0 + f1;
  j[5]  = j0 + f2;
  j[6]  = j0 + f3;
  j[7]  = j0 + f4;
  j[8]  = j0 + f1;
  j[9]  = j0 + f2;
  j[10] = j0 + f3;
  j[11] = j0 + f4;
  j[12] = j0 + f1;
  j[13] = j0 + f2;
  j[14] = j0 + f3;
  j[15] = j0 + f4;
  if( dy<0 ){
    f1 = -1;
    f2 = 0;
    f3 = 1;
    f4 = 2;
  } else {
    f1 = 2;
    f2 = 1;
    f3 = 0;
    f4 = -1;
  }
  i[0]  = i0 + f1;
  i[1]  = i0 + f1;
  i[2]  = i0 + f1;
  i[3]  = i0 + f1;  
  i[4]  = i0 + f2;
  i[5]  = i0 + f2;
  i[6]  = i0 + f2;
  i[7]  = i0 + f2;  
  i[8]  = i0 + f3;
  i[9]  = i0 + f3;
  i[10] = i0 + f3;
  i[11] = i0 + f3;  
  i[12] = i0 + f4;
  i[13] = i0 + f4;
  i[14] = i0 + f4;
  i[15] = i0 + f4;  
  return true;
}


double RectGrid::interp2d_nearest(double x,double y,double* f0){
  if( !this->point_between_pixel_centers(x,y,0) ){
    return 0.0;
  }

  int i0,j0;
  this->match_point_to_pixel(x,y,i0,j0);
  return f0[i0*this->Nx+j0];
}

double RectGrid::interp2d_bilinear(double x,double y,double* f0){
   if( !this->point_between_pixel_centers(x,y,0) ){
     return 0.0;
   }

  //    a     b
  // 0-----------1
  // | w3  |  w2 |  c
  // |---(x,y)---|
  // | w1  |  w0 |  d
  // 2-----------3
  //
  // w3=ac, w2=bc, w1=ad, w0=bd
  
  int ii[4],jj[4];
  match_point_to_closest_4(x,y,ii,jj);

  double a = (x - this->center_x[jj[0]]);
  double b = (this->center_x[jj[1]] - x);
  double c = (this->center_y[ii[0]] - y);
  double d = (y - this->center_y[ii[2]]);

  std::vector<double> w{b*d,a*d,b*c,a*c};

  double z_interp = 0.0;
  for(int k=0;k<4;k++){
    z_interp += w[k]*f0[ii[k]*this->Nx+jj[k]];
  }
  return z_interp/(this->step_x*this->step_y);
}

double RectGrid::interp2d_bicubic(double x,double y,double* dummy){
  if( !this->point_between_pixel_centers(x,y,1) ){
    return 0.0;
  }
  
  if( this->zx == NULL ){
    this->calculate_zx();
  }
  if( this->zy == NULL ){
    this->calculate_zy();
  }
  if( this->zxy == NULL ){
    this->calculate_zxy();
  }

  int ii[4],jj[4];
  match_point_to_closest_4(x,y,ii,jj);

  double tab_left[16]  = {1,0,0,0,0,0,1,0,-3,3,-2,-1,2,-2,1,1};
  double tab_right[16] = {1,0,-3,2,0,0,3,-2,0,1,-2,1,0,0,-1,1};

  int x00 = ii[2]*this->Nx+jj[2];
  int x01 = ii[0]*this->Nx+jj[0];
  int x10 = ii[3]*this->Nx+jj[3];
  int x11 = ii[1]*this->Nx+jj[1];
  double tab_mid[16] = {this->z[x00],this->z[x01],this->step_y*this->zy[x00],this->step_y*this->zy[x01],this->z[x10],this->z[x11],this->step_y*this->zy[x10],this->step_y*this->zy[x11],this->step_x*this->zx[x00],this->step_x*this->zx[x01],this->step_y*this->step_x*this->zxy[x00],this->step_y*this->step_x*this->zxy[x01],this->step_x*this->zx[x10],this->step_x*this->zx[x11],this->step_y*this->step_x*this->zxy[x10],this->step_y*this->step_x*this->zxy[x11]};
  //  double tab_mid[16] = {this->z[x00],this->z[x01],this->zy[x00],this->zy[x01],this->z[x10],this->z[x11],this->zy[x10],this->zy[x11],this->zx[x00],this->zx[x01],this->zxy[x00],this->zxy[x01],this->zx[x10],this->zx[x11],this->zxy[x10],this->zxy[x11]};

  double coeff[16],tab_tmp[16];
  this->multiply_table_table(4,4,4,tab_mid,tab_right,tab_tmp);
  this->multiply_table_table(4,4,4,tab_left,tab_tmp,coeff);

  double xu = (x-this->center_x[jj[0]])/this->step_x;
  double yu = (y-this->center_y[ii[2]])/this->step_y;

  double vec_left[4]  = {1,xu,pow(xu,2),pow(xu,3)};
  double vec_right[4] = {1,yu,pow(yu,2),pow(yu,3)};
  double vec_tmp[4];
  this->multiply_table_vector(4,4,coeff,vec_right,vec_tmp);
  return multiply_vector_vector(4,vec_left,vec_tmp);
}


void RectGrid::calculate_zx(){
  this->zx = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_1(this->Nx,this->Ny,this->center_x,this->z,this->zx);
}

void RectGrid::calculate_zy(){
  this->zy = (double*) calloc(this->Nz,sizeof(double));
  double* tmp_1 = (double*) calloc(this->Nz,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zy[j*this->Ny+i] = this->z[i*this->Nx+j];
    }
  }
  this->calculate_derivative_1(this->Ny,this->Nx,this->center_y,this->zy,tmp_1);
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zy[i*this->Nx+j] = tmp_1[j*this->Ny+i];
    }
  }
  free(tmp_1);
}

void RectGrid::calculate_zxy(){
  if( this->zy == NULL ){
    this->calculate_zy();
  }
  this->zxy = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_1(this->Nx,this->Ny,this->center_x,this->zy,this->zxy);  
}

void RectGrid::calculate_zxx(){
  this->zxx = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_2(this->Nx,this->Ny,this->center_x,this->z,this->zxx);
}

void RectGrid::calculate_zyy(){
  this->zyy = (double*) calloc(this->Nz,sizeof(double));
  double* tmp_1 = (double*) calloc(this->Nz,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zyy[j*this->Ny+i] = this->z[i*this->Nx+j];
    }
  }
  this->calculate_derivative_2(this->Ny,this->Nx,this->center_y,this->zyy,tmp_1);
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zyy[i*this->Nx+j] = tmp_1[j*this->Ny+i];
    }
  }
  free(tmp_1);
}


void RectGrid::calculate_derivative_1(int Nh,int Nv,double* h,double* zz,double* zout){
  // calculate the FIRST derivative along the horizontal axis
  int h0,v0;
  std::vector<int> rel_index_v;
  std::vector<int> rel_index_h;
  std::vector<double> coeff;
  double dh = h[1]-h[0];
  
  // 1st column
  if( this->options["dev1_accu"] == "1" ){
    rel_index_v = {0,0};
    rel_index_h = Constants::derivative_1_forward_1_index;
    coeff       = Constants::derivative_1_forward_1_coeff;
  } else {
    rel_index_v = {0,0,0};
    rel_index_h = Constants::derivative_1_forward_2_index;
    coeff       = Constants::derivative_1_forward_2_coeff;
  }
  h0 = 0;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
  }

  // last column
  if( this->options["dev1_accu"] == "1" ){
    rel_index_v = {0,0};
    rel_index_h = Constants::derivative_1_backward_1_index;
    coeff       = Constants::derivative_1_backward_1_coeff;
  } else {
    rel_index_v = {0,0,0};
    rel_index_h = Constants::derivative_1_backward_2_index;
    coeff       = Constants::derivative_1_backward_2_coeff;
  }
  h0 = Nh-1;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
  }

  // middle chunk
  rel_index_v = {0,0,0};
  rel_index_h = Constants::derivative_1_central_2_index;
  coeff       = Constants::derivative_1_central_2_coeff;
  for(int v=0;v<Nv;v++){
    for(int h=1;h<Nh-1;h++){
      zout[v*Nh+h] = this->weighted_sum(v,h,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
    }
  }  
}


void RectGrid::calculate_derivative_2(int Nh,int Nv,double* h,double* zz,double* zout){
  // calculate the SECOND derivative along the horizontal axis
  int h0,v0;
  std::vector<int> rel_index_v;
  std::vector<int> rel_index_h;
  std::vector<double> coeff;
  double dh2 = pow(h[1]-h[0],2);
  
  // 1st column
  if( this->options["dev2_accu"] == "1" ){
    rel_index_v = {0,0,0};
    rel_index_h = Constants::derivative_2_forward_1_index;
    coeff       = Constants::derivative_2_forward_1_coeff;
  } else {
    rel_index_v = {0,0,0,0};
    rel_index_h = Constants::derivative_2_forward_2_index;
    coeff       = Constants::derivative_2_forward_2_coeff;
  }
  h0 = 0;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
  }

  // last column
  if( this->options["dev2_accu"] == "1" ){
    rel_index_v = {0,0,0};
    rel_index_h = Constants::derivative_2_backward_1_index;
    coeff       = Constants::derivative_2_backward_1_coeff;
  } else {
    rel_index_v = {0,0,0,0};
    rel_index_h = Constants::derivative_2_backward_2_index;
    coeff       = Constants::derivative_2_backward_2_coeff;
  }
  h0 = Nh-1;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
  }

  // middle chunk
  rel_index_v = {0,0,0};
  rel_index_h = Constants::derivative_2_central_2_index;
  coeff       = Constants::derivative_2_central_2_coeff;
  for(int v=0;v<Nv;v++){
    for(int h=1;h<Nh-1;h++){
      zout[v*Nh+h] = this->weighted_sum(v,h,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
    }
  }  
}

double RectGrid::weighted_sum(int i0,int j0,std::vector<int> rel_i,std::vector<int> rel_j,std::vector<double> coeff,int z_Nx,double* zz){
  double sum = 0.0;
  for(int k=0;k<rel_i.size();k++){
    int index_z  = (i0+rel_i[k])*z_Nx + (j0+rel_j[k]);
    sum += coeff[k]*zz[index_z];
  }
  return sum;
}

void RectGrid::multiply_vector_scalar(std::vector<double> &v,double k){
  std::transform(v.begin(),v.end(),v.begin(), [k](double &c){ return c*k; });
}

void RectGrid::multiply_table_vector(int Nrows,int Ncols,double* tab_input,double* vec_input,double* vec_output){
  for(int i=0;i<Nrows;i++){
    double sum = 0.0;
    for(int j=0;j<Ncols;j++){
      sum += tab_input[i*Nrows+j] * vec_input[j];
    }
    vec_output[i] = sum;
  }
}

void RectGrid::multiply_table_table(int Nrows_1,int Ncols_1,int Ncols_2,double* tab_input_1,double* tab_input_2,double* tab_output){
  for(int i=0;i<Nrows_1;i++){
    for(int k=0;k<Ncols_2;k++){
      double sum = 0.0;
      for(int j=0;j<Ncols_1;j++){
	sum += tab_input_1[i*Ncols_1+j] * tab_input_2[j*Ncols_2+k];
      }
      tab_output[i*Ncols_2+k] = sum;
    }
  }
}
  
double RectGrid::multiply_vector_vector(int N,double* vec_1,double* vec_2){
  double sum = 0.0;
  for(int i=0;i<N;i++){
    sum += vec_1[i]*vec_2[i];
  }
  return sum;
}

