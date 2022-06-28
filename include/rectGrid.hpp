#ifndef RECT_GRID_HPP
#define RECT_GRID_HPP

#include <string>
#include <vector>
#include <map>

namespace gerlumph {

  class Constants{
  public:
    // Finite difference coefficients and relative indices.
    // 1st derivatives
    static const std::vector<int>    derivative_1_forward_1_index;
    static const std::vector<double> derivative_1_forward_1_coeff;
    static const std::vector<int>    derivative_1_forward_2_index;
    static const std::vector<double> derivative_1_forward_2_coeff;
    static const std::vector<int>    derivative_1_backward_1_index;
    static const std::vector<double> derivative_1_backward_1_coeff;
    static const std::vector<int>    derivative_1_backward_2_index;
    static const std::vector<double> derivative_1_backward_2_coeff;
    static const std::vector<int>    derivative_1_central_2_index;
    static const std::vector<double> derivative_1_central_2_coeff;
    // 2nd derivatives
    static const std::vector<int>    derivative_2_forward_1_index;
    static const std::vector<double> derivative_2_forward_1_coeff;
    static const std::vector<int>    derivative_2_forward_2_index;
    static const std::vector<double> derivative_2_forward_2_coeff;
    static const std::vector<int>    derivative_2_backward_1_index;
    static const std::vector<double> derivative_2_backward_1_coeff;
    static const std::vector<int>    derivative_2_backward_2_index;
    static const std::vector<double> derivative_2_backward_2_coeff;
    static const std::vector<int>    derivative_2_central_2_index;
    static const std::vector<double> derivative_2_central_2_coeff;
  };


  class RectGrid{
  public:
    int Nx;
    int Ny;
    int Nz;              // Total number of pixels in the grid
    double* center_x;    // 1D, size of Nx: pixel center x coordinates
    double* center_y;    // 1D, size of Ny: pixel center y coordinates
    double* z   = NULL;  // 1D, size of Nx*Ny: Values of the 2D surface in row-major format, i.e. z[i*Nx+j]
    double* zx  = NULL;  // same as z: first derivative along the x axis
    double* zy  = NULL;  // same as z: first derivative along the y axis
    double* zxy = NULL;  // same as z: mixed derivative along the x and y axes
    double* zxx = NULL;  // same as z: second derivative along the x axis
    double* zyy = NULL;  // same as z: second derivative along the y axis
  
    double* bound_x; // 1D, size of Nx+1: pixel x boundaries, can be accessed by x[j],x[j+1], where 0<j<Nx
    double* bound_y; // 1D, size of Ny+1: pixel y boundaries, can be accessed by y[i],y[i+1], where 0<i<Ny
    double width;
    double height;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double step_x;
    double step_y;

    std::map<std::string,std::string> options;
  
    RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::map<std::string,std::string> options = {{"dev1_accu","1"},{"dev2_accu","2"},{"interp","bilinear"}});
    RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,const std::string filepath,std::map<std::string,std::string> options = {{"dev1_accu","1"},{"dev2_accu","2"},{"interp","bilinear"}});
    RectGrid(const RectGrid& grid);
    ~RectGrid();

    void set_interp(std::string interp);
    RectGrid embeddedNewGrid(int new_Nx,int new_Ny,std::string mode="interp");
  
    bool point_in_grid(double x,double y);
    bool point_between_pixel_centers(double x,double y,int boundary_size);
    void match_point_to_pixel(double x,double y,int& i0,int& j0);
    bool match_point_to_closest_4(double x,double y,int* i,int* j);
    bool match_point_to_closest_16(double x,double y,int* i,int* j);

    double (RectGrid::*interp2d)(double x,double y,double* f0);
    double interp2d_nearest(double x,double y,double* f0);
    double interp2d_bilinear(double x,double y,double* f0);
    double interp2d_bicubic(double x,double y,double* dummy);
  
    void calculate_zx();
    void calculate_zy();
    void calculate_zxy();
    void calculate_zxx();
    void calculate_zyy();

  private:
    void common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::map<std::string,std::string> options);
  
    void multiply_vector_scalar(std::vector<double> &v,double k);
    void multiply_table_vector(int Nrows,int Ncols,double* tab_input,double* vec_input,double* vec_output);
    void multiply_table_table(int Nrows_1,int Ncols_1,int N_cols_2,double* tab_input_1,double* tab_input_2,double* tab_output);
    double multiply_vector_vector(int N,double* vec_1,double* vec_2);

    void calculate_derivative_1(int Nh,int Nv,double* h,double* zz,double* zout);
    void calculate_derivative_2(int Nh,int Nv,double* h,double* zz,double* zout);

    double weighted_sum(int i0,int j0,std::vector<int> rel_i,std::vector<int> rel_j,std::vector<double> coeff,int z_Nx,double* zz);
  };

}
  
#endif /* RECT_GRID_HPP */
