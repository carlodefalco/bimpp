/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>
 
#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>

#include "Taylor_Disc_Galerkin.h"


// mpirun -np 4 mainDG $PWD inputs/dem_real_vladi.octbin.gz inputs/mask_in_vladi.octbin.gz 
// mpirun -np 4 mainDG $PWD inputs/dem_ideal.octbin.gz inputs/mask_in_vladi.octbin.gz 

static constexpr char VARNAME_1[255] = "dem";
static constexpr char VARNAME_2[255] = "mask_in";

// properties of the input dem
static constexpr double res = 2e-2;//0.005*500; // it is also the minimum resolution of the bim element
static constexpr double Nx = 101;//101;//165;//201;//188; // # columns
static constexpr double Ny = 101;//101;//175;//201;//180; // # rows

 
static constexpr double L = res*(Nx-1);
static constexpr double H = res*(Ny-1);
static std::vector<double>   dem;
static std::vector<double>   dem_slope_x;
static std::vector<double>   dem_slope_y;
static std::vector<double>   basin_mask;
static constexpr int NUM_REFINEMENTS  = 6; // 8
static constexpr int NUM_TREFINEMENTS = 1; // 10



static constexpr double SPACE_ADAPTDT = 4e-2;//1e-2; // put zero if you want at each time step
static constexpr double SAVEDT = 1.e-2; // must never be null 
static constexpr double DELTAT = 1.e-3;
static constexpr double REDCDT = .5; 
static constexpr double T      = .5; 
 
static constexpr bool is_time_adaptivity    = false; 
static constexpr bool is_initial_refinement = false;
static constexpr bool is_space_adaptivity   = false;
static constexpr bool is_non_reflBC         = true; 
static constexpr bool is_bed_friction       = false; 
static constexpr bool is_stress_tensor      = false;



static constexpr double h_min = 1e-5;
static constexpr double grav = 9.81;
static constexpr double density = 1291.;
static constexpr double turbulence_coeff = 1.e8;
static constexpr double surface_pressure = 0;//101325.;
static constexpr double bed_friction_angle_rad = 22.*M_PI/180; //33.9*M_PI/180; //0.0; //23*M_PI/180; 
static constexpr double fluid_viscosity = 5e1;
static constexpr double yield_shear_stress = 2e3;//.5*density*grav*38*std::sin(bed_friction_angle_rad);

static constexpr double level_wet           = 3;  
static constexpr double level_interface     = 6; // minimum resolution! 
static constexpr double mesh_size_dry       = res*100;//res/60*std::pow(2,level_interface); //res*std::pow(2,level_interface); 
static constexpr double mesh_size_wet       = res/2;///10;//res/20;//res;//mesh_size_dry/std::pow(2,level_wet); // finest resolution
static constexpr double mesh_size_interface = res/20;//res/30;//res/60;//mesh_size_dry/std::pow(2,level_interface);
 

// Connectivity of local element
constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0,  0,
   0,  H,
   L,  0,
   L,  H};

const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {  1,    3,    4,    2,    1 };

// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
hanging_refinement (tmesh::quadrant_iterator quadrant)
{
  
  double x_minus, x_plus, y_minus, y_plus;
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  double x_center, y_center;
  x_center = quadrant->centroid (0);
  y_center = quadrant->centroid (1);
  
  const auto marker = x_center<3./4.*L && x_center>L/4. && y_center<3./4.*H && y_center>H/4.;
//  const auto marker = (x_minus+x_plus)/2<L/2;//(x_minus+x_plus)/2<L/2 && (y_minus+y_plus)/2>H/2 ? 1 : 0;
  return marker; 
}


static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}

static std::array<int,3>
global_coord_2_raster(const double& x,
                      const double& y)
{
  static double i_x;
  static double i_y;
  static int ii;
  
  
  // nearest neighbor
  i_x =   std::round(x / res);
  i_y = - std::round(y / res) + (Ny-1);
  
  ii = raster_2_vector(i_x,i_y);
  
  return(std::array<int,3>{{ ii,int(i_x),int(i_y) }});
}



using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>; //distributed_vector; //std::vector<double>;         // Typedef for local q_0 vector // distributed_vector

//double h0_fun (const double& xx, const double& yy)  { return std::max (0., (8. - std::sin (M_PI * xx / 2. / 400.) - dem[global_coord_2_raster(xx,yy)[0]])); }
double h0_fun (const double& xx, const double& yy) 
{
  return ( 1.+1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ) );
  //return( std::abs(xx-L/2.)<=L/2. && std::abs(yy-H/2.)<=L.2 ? 2 : 1. ); 
  //return(std::sqrt(std::pow(xx-L/2.,2.) + std::pow(yy-H/2.,2.))<=.25*std::min(L,H) ? 2. : 1. ); 
  return(xx<=L/2. ? 2. : 1. ); 


  // const double HH = 1.;
  // const double omega = (std::pow((xx-.5*L)/L,2.) + std::pow((yy-.5*L)/L,2.)) <= std::pow((.2 + .01 * std::sin(10.*M_PI*(yy-.5*L)/L)),2.) ? 1. : 0.;
  // return(std::max (0., std::min (60.-(100 - 100 * xx/L), HH)) * omega); 
  
  
  //return(std::sqrt(std::pow(xx-L/2.,2.) + std::pow(yy-H/2.,2.))<=150 ? 70 : 0. ); 

  return(basin_mask[global_coord_2_raster(xx,yy)[0]]==1 ? 38 : 0.); 
  //return (xx<=L/2. ? 70 : 7.); //(xx<=L/2. ? 70 : 0.);
  return ( 1.+1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ) );
  //return ( 0.+1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ) );

  
  
  if (xx>L*1./4. && xx<L*3./4. && yy >H*1./4. && yy <H*3./4.) //(xx>L*3./10. && xx<7./10.*L)
  {
    //return (1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ));
    return 40.;
  } 
  
//  if (xx>L/4 && xx<3/4*L && yy>H/4 && yy <3/4*H)
//  {
//    return 2;
//  }
  return 0.;
} 
double Ux0_fun (double xx, double yy) { return 0.; }
double Uy0_fun (double xx, double yy) { return 0.; }


// Assemble vector from mesh.
// FIXME  the following two functions are copied over from
// "quad_operators.cpp" as they were not exported in an header,
// should find better way to avoid code duplication
void
assemble_vector (tmesh::quadrant_iterator& quadrant,
                 const std::array<double, 4>& locrhs,
                 Q1& rhs,
                 const ordering& ord = default_ord)
{
  
  std::vector<unsigned int> rows;
  rows.reserve (2);
  int i, r;
  
  for (i = 0; i < 4; ++i)
  {
    rows.clear ();
    
    if (! quadrant->is_hanging (i))
      rows.push_back (quadrant->gt (i));
    else
    {
      rows.push_back (quadrant->gparent (0, i));
      rows.push_back (quadrant->gparent (1, i));
    }
    
    for (r = 0; r < rows.size (); ++r)
    rhs[ord (rows[r])] +=
    locrhs[i] / rows.size ();
  }
}


void
compute_slope()
{

  int ii, jj;
  for (ii=1; ii<(Ny-1); ii++)
  {
    for (jj=1; jj<(Nx-1); jj++)
    {
      const auto i_vec = raster_2_vector(jj,ii);

      const auto i_vec_north = raster_2_vector(jj,ii-1);
      const auto i_vec_south = raster_2_vector(jj,ii+1);

      const auto i_vec_east = raster_2_vector(jj+1,ii);
      const auto i_vec_west = raster_2_vector(jj-1,ii);

      dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec_west ])/(2*res);
      dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec_south])/(2*res);
    }
  }

  ii = 0;
  for (jj=1; jj<(Nx-1); jj++)
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_south = raster_2_vector(jj,ii+1);

    const auto i_vec_east = raster_2_vector(jj+1,ii);
    const auto i_vec_west = raster_2_vector(jj-1,ii);

    dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec_west ])/(2*res);
    dem_slope_y[i_vec] = (dem[i_vec]-dem[i_vec_south])/res;
  }

  ii = Ny-1;
  for (jj=1; jj<(Nx-1); jj++)
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_north = raster_2_vector(jj,ii-1);

    const auto i_vec_east = raster_2_vector(jj+1,ii);
    const auto i_vec_west = raster_2_vector(jj-1,ii);

    dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec_west ])/(2*res);
    dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec])/res;
  }

  jj = 0;
  for (ii=1; ii<(Ny-1); ii++)
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_north = raster_2_vector(jj,ii-1);
    const auto i_vec_south = raster_2_vector(jj,ii+1);

    const auto i_vec_east = raster_2_vector(jj+1,ii);

    dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec ])/res;
    dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec_south])/(2*res);
  }

  jj = Nx-1;
  for (ii=1; ii<(Ny-1); ii++)
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_north = raster_2_vector(jj,ii-1);
    const auto i_vec_south = raster_2_vector(jj,ii+1);

    const auto i_vec_west = raster_2_vector(jj-1,ii);

    dem_slope_x[i_vec] = (dem[i_vec ]-dem[i_vec_west ])/res;
    dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec_south])/(2*res);
  }


  // compute corner points
  ii = 0; jj = 0;
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_south = raster_2_vector(jj,ii+1);

    const auto i_vec_east = raster_2_vector(jj+1,ii);

    dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec])/res;
    dem_slope_x[i_vec] = (dem[i_vec]-dem[i_vec_south])/res;
  }

  ii = Ny-1; jj = 0;
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_north = raster_2_vector(jj,ii-1);

    const auto i_vec_east = raster_2_vector(jj+1,ii);

    dem_slope_x[i_vec] = (dem[i_vec_east ]-dem[i_vec])/res;
    dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec])/res;
  }

  ii = 0; jj = Nx-1;
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_south = raster_2_vector(jj,ii+1);

    const auto i_vec_west = raster_2_vector(jj-1,ii);

    dem_slope_x[i_vec] = (dem[i_vec]-dem[i_vec_west ])/res;
    dem_slope_y[i_vec] = (dem[i_vec]-dem[i_vec_south])/res;
  }  

  ii = Ny-1; jj = Nx-1;
  {
    const auto i_vec = raster_2_vector(jj,ii);

    const auto i_vec_north = raster_2_vector(jj,ii-1);

    const auto i_vec_west = raster_2_vector(jj-1,ii);

    dem_slope_x[i_vec] = (dem[i_vec]-dem[i_vec_west ])/res;
    dem_slope_y[i_vec] = (dem[i_vec_north]-dem[i_vec])/res;
  }

}



template <class T>
void
quadrant_marker_list (tmesh::quadrant_iterator& q, 
                      const T& only_h, const T& only_Ux, const T& only_Uy, const double& dt, std::set<int>& output)
{

  std::array<double,4> h_current = {0,0,0,0};
  std::array<double,2> vel       = {0,0};
  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii)){
        const auto & h_candidate = only_h[q->gt (ii)];
        h_current[ii] = h_candidate;

        vel[0] += h_candidate>h_min ? only_Ux[q->gt (ii)]/h_candidate : 0.;
        vel[1] += h_candidate>h_min ? only_Uy[q->gt (ii)]/h_candidate : 0.;
      }
      else
      {
        const auto & h_candidate = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        h_current[ii] = h_candidate;

        vel[0] += h_candidate>h_min ? .5*(only_Ux[q->gparent (0,ii)] + only_Ux[q->gparent (1,ii)])/h_candidate : 0.;

        vel[1] += h_candidate>h_min ? .5*(only_Uy[q->gparent (0,ii)] + only_Uy[q->gparent (1,ii)])/h_candidate : 0.;
      }
    }
  vel[0] /= 4.;
  vel[1] /= 4.;

  const bool basin_check = ((h_current[0]+h_current[1]+h_current[2]+h_current[3])> h_min && 
                            (h_current[0]*h_current[1]*h_current[2]*h_current[3])<=h_min) ? true : false;


  if (basin_check)
  {
    std::vector<std::tuple<tmesh::quadrant_iterator, int> > quadrant_list;
    output.insert(q->get_global_quad_idx ());
    for (auto quadrant_nei  = q->begin_neighbor_sweep();
              quadrant_nei != q->end_neighbor_sweep (); ++quadrant_nei)
    {
      quadrant_list.push_back(std::tuple<tmesh::quadrant_iterator, int>{quadrant_nei, quadrant_nei->get_global_quad_idx ()});
    }

    const double xx_ini = q->centroid (0);
    const double yy_ini = q->centroid (1);

    const double xx_fin = xx_ini + vel[0]*dt;
    const double yy_fin = yy_ini + vel[1]*dt;

    for (auto & qq : quadrant_list)
    {
      auto & quadrant = std::get<0>(qq);

      //std::cout << quadrant->p(0,0) << " " << std::get<1>(qq) << std::endl;

      const double x1 = quadrant->p(0,0);
      const double x2 = quadrant->p(0,1);
      const double y1 = quadrant->p(1,0);
      const double y2 = quadrant->p(1,2); 

      //std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;

      const bool cond1 = ( yy_fin - y1)>=0;
      const bool cond2 = (-yy_fin + y2)>=0;
      const bool cond3 = ( xx_fin - x1)>=0;
      const bool cond4 = (-xx_fin + x2)>=0;


      if (cond1 && cond2 && cond3 && cond4) // the point is internal to the current quadrant
      {
        output.insert(std::get<1>(qq));
      }

    }
  }

}







// Re-Define tic and toc to add an MPI_Barrier
#define TIC()  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int
main (int argc, char **argv)
{
  // Management of solutions ordering
  ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<12, 0> (gt); };
  ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<12, 4> (gt); };
  ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<12, 8> (gt); };

  // Management of solutions ordering
  ordering ordh_nodal  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 0> (gt); };
  ordering ordUx_nodal = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 1> (gt); };
  ordering ordUy_nodal = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 2> (gt); };
  
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (argc != 4) 
  {
    std::cerr << "You should provide as input respectively, save directory, dem directory, mask directory" << std::endl;

    // Close MPI and print report
    MPI_Barrier (MPI_COMM_WORLD);
    if (rank == 0) { print_timing_report (); }
    MPI_Finalize ();
    return 0;
  }

  const auto SAVE_DIR = argv[1];
  const auto DEM_DIR  = argv[2]; 
  const auto MASK_DIR = argv[3];

  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  
  
  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  
  // for (int ii=0; ii<1; ii++)
  // {
  //   tmsh.set_refine_marker (hanging_refinement);
  //   tmsh.refine (recursive, 1);
  // }
  TOC ("Uniform refinement");
  
  
  
  // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes (); // Return total number of nodes owned by all process
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes (); // Return number of nodes owned by local process
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();  // Return number of quadrants owned by local process across all trees
  tmesh::idx_t gn_elements = tmsh.num_global_quadrants (); // Return number of quadrants owned by all processes across all trees
  

  /// Allocate initial data container

  auto elemental_mass_matrix_inverted = [] (tmesh::quadrant_iterator quadrant)
  {

/*

close all
clear
%clc

wGP = .25*ones(4,1);
xiGP = [(-1/sqrt(3)+1)/2, (-1/sqrt(3)+1)/2, (1/sqrt(3)+1)/2, (1/sqrt(3)+1)/2];
yiGP = [(-1/sqrt(3)+1)/2, (1/sqrt(3)+1)/2, (-1/sqrt(3)+1)/2, (1/sqrt(3)+1)/2];

figure()
plot(xiGP, yiGP, 'o')
axis([0 1 0 1])


Me = zeros(4,4);
for j=1:4
    phi = BaseFunc(xiGP(j), yiGP(j));
    for k=1:4
        for l=1:4
            Me(k,l) = Me(k,l) +wGP(j)*phi(k)*phi(l);
        end
    end 
end
% Inverse
iMe = inv(Me)


%% utils

function phi = BaseFunc(x, y)
% order as bimpp  
phi(4) = x*y;
phi(3) = (1-x)*y;
phi(2) = (1-y)*x;
phi(1) = (1-x)*(1-y);
end



*/




    std::array<double,4> isdof_or_hanging;

    double Dx, Dy, area;
    Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
    Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
    area = Dx * Dy;

    for (int ii = 0; ii < 4; ++ii){
      if (! quadrant->is_hanging (ii)){
        isdof_or_hanging[ii] = 1.;
      } else {
        isdof_or_hanging[ii] = .5;
      }
    }

    std::array<std::array<double,4>, 4> iMe = {std::array<double,4>{{16., -8., -8., 4.}}, 
                                               std::array<double,4>{{-8., 16., 4., -8.}},  
                                               std::array<double,4>{{-8., 4., 16., -8.}}, 
                                               std::array<double,4>{{4., -8., -8., 16.}}};

    for (int ii = 0; ii < 4; ++ii)
    {
      for (int jj = 0; jj < 4; ++jj)
      {
        iMe[ii][jj] /= area*isdof_or_hanging[ii]*isdof_or_hanging[jj];
      }
    }

    return (iMe); 
  };

  Q0 sol_onehalf (ln_elements * 3);
  sol_onehalf.assign(sol_onehalf.size(), 0.0);

  distributed_vector dg_coefficients (ln_elements * 4 * 3), incr_dg_coefficients(ln_elements * 4 * 3);
  dg_coefficients.get_owned_data ().assign(dg_coefficients.size(), 0.);
  incr_dg_coefficients.get_owned_data ().assign(incr_dg_coefficients.size(), 0.);
  
  // for visualization
  Q1 sol  (ln_nodes * 3);
  sol.get_owned_data  ().assign (sol.get_owned_data  ().size (), 0.0);

  Q1 soll(ln_nodes * 3);
  soll.get_owned_data ().assign(soll.size(), 0.);

  Q1 Z (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);

  Q0 slope_x (ln_elements);
  slope_x.assign (slope_x.size (), 0.0);

  Q0 slope_y (ln_elements);
  slope_y.assign (slope_y.size (), 0.0);

  std::string str = ""; 
  char filename[255]="", arr[255]="";

  TIC();
  str = std::string(DEM_DIR); 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);


  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
  octave_value v;

  octave_io_open (filename, m_in, &m_out);
  octave_load (VARNAME_1, v);
  Matrix M = v.matrix_value ();
  dem.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());

  str = std::string(MASK_DIR); 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);

  octave_io_open (filename, m_in, &m_out);
  octave_load (VARNAME_2, v);
  M = v.matrix_value ();
  basin_mask.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), basin_mask.begin ());
  TOC("Load data matrix");
  
  // compute raster slope
  dem_slope_x.resize(Nx*Ny);
  dem_slope_y.resize(Nx*Ny);
  compute_slope();



  // Initialize 
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
  {

    const auto iMe = elemental_mass_matrix_inverted(quadrant);

    const auto Dx = quadrant->p(0, 1) - quadrant->p(0, 0);
    const auto Dy = quadrant->p(1, 2) - quadrant->p(1, 0);
    const auto area = Dx * Dy;

    double xx_c=quadrant->centroid(0);
    double yy_c=quadrant->centroid(1); 

    slope_x[quadrant->get_forest_quad_idx ()] = dem_slope_x[global_coord_2_raster(xx_c,yy_c)[0]];
    slope_y[quadrant->get_forest_quad_idx ()] = dem_slope_y[global_coord_2_raster(xx_c,yy_c)[0]];

    std::array<double,4> dg_coeff_h = {0., 0., 0., 0.}, dg_coeff_Ux = {0., 0., 0., 0.}, dg_coeff_Uy = {0., 0., 0., 0.};
    for (int kk = 0; kk < 4; ++kk)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii); 

        // trapezoidal rule
        if (! quadrant->is_hanging (ii)){
          dg_coeff_h [kk] += .25 * area * h0_fun   (xx, yy) * 1. * (ii == kk);
          dg_coeff_Ux[kk] += .25 * area * Ux0_fun  (xx, yy) * 1. * (ii == kk);
          dg_coeff_Uy[kk] += .25 * area * Uy0_fun  (xx, yy) * 1. * (ii == kk);
        }
        else
        {
          dg_coeff_h [kk] += .25 * area * h0_fun   (xx, yy) * .5 * (ii == kk);
          dg_coeff_Ux[kk] += .25 * area * Ux0_fun  (xx, yy) * .5 * (ii == kk);
          dg_coeff_Uy[kk] += .25 * area * Uy0_fun  (xx, yy) * .5 * (ii == kk);
        }
        
      }
    }

    std::array<double,4> dg_coeff_h_, dg_coeff_Ux_, dg_coeff_Uy_;
    for (int kk = 0; kk < 4; ++kk)
    {
      const auto & iMe_kk = iMe[kk];

      const double sum_iMe_kk = std::accumulate(iMe_kk.begin(), iMe_kk.end(), 0.);

      // matrice consistente
      // dg_coeff_h_ [kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_h .begin(), 0.); 
      // dg_coeff_Ux_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_Ux.begin(), 0.); 
      // dg_coeff_Uy_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_Uy.begin(), 0.); 

      //std::cout << iMe[kk][0] << " " << iMe[kk][1] << " " << iMe[kk][2] << " " << iMe[kk][3] << " " << sum_iMe_kk << std::endl;

      // matrice lumpata
      dg_coeff_h_ [kk] = sum_iMe_kk*dg_coeff_h [kk]; 
      dg_coeff_Ux_[kk] = sum_iMe_kk*dg_coeff_Ux[kk]; 
      dg_coeff_Uy_[kk] = sum_iMe_kk*dg_coeff_Uy[kk]; 
    }

    for (int kk = 0; kk < 4; ++kk)
    {
      // if (dg_coeff_h_ [kk]!=0)
      // std::cout << dg_coeff_h_ [0] << " " << dg_coeff_h_ [1] << " " << dg_coeff_h_ [2] << " " << dg_coeff_h_ [3] << std::endl;

      dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_h_ [kk];
      dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_Ux_[kk];
      dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_Uy_[kk];

      incr_dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + kk] = 0.;
      incr_dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + kk] = 0.;
      incr_dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + kk] = 0.;
    }
    


    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {
      for (int kk = 0; kk < 4; ++kk)
      {
        dg_coefficients[ordh  (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
        dg_coefficients[ordUx (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
        dg_coefficients[ordUy (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;

        incr_dg_coefficients[ordh  (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
        incr_dg_coefficients[ordUx (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
        incr_dg_coefficients[ordUy (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
      }
    }



    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii); 
        
        Z[quadrant->gt (ii)] = dem[global_coord_2_raster(xx,yy)[0]]; 
      }
      
      else
      {
        // touch parent nodes to set up distributed vector structure
        Z   [quadrant->gparent(0,ii)] += 0.;
        Z   [quadrant->gparent(1,ii)] += 0.;
      }
    }
  }


  // bim2a_solution_with_ghosts in quad_operators.cpp
  dg_coefficients.remap();
  dg_coefficients.assemble(replace_op);

  Z.remap();
  Z.assemble(replace_op);

  incr_dg_coefficients.remap();
  incr_dg_coefficients.assemble(replace_op);


  TIC();
  Q1 mass (ln_nodes * 3);
  bim2a_mass_vector (tmsh, mass, ordh_nodal);
  bim2a_mass_vector (tmsh, mass, ordUx_nodal);
  bim2a_mass_vector (tmsh, mass, ordUy_nodal);
  mass.assemble ();


  for (auto quadrant = tmsh.begin_quadrant_sweep ();
   quadrant != tmsh.end_quadrant_sweep ();
   ++quadrant)
  {

    std::array<double,4> isdof_or_hanging;

    const auto Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
    const auto Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
    const auto area = Dx * Dy;


    for (int ii = 0; ii < 4; ++ii){
      if (! quadrant->is_hanging (ii)){
        isdof_or_hanging[ii] = 1.;
      } else {
        isdof_or_hanging[ii] = .5;
      }
    }

    // midpoint integration
    double inh = -1.e4;
    double integral_h = 0., integral_Ux = 0., integral_Uy = 0.;
    for (int ii = 0; ii < 4; ++ii)
    {
      inh  = std::max(inh,  dg_coefficients[ordh    (quadrant->get_global_quad_idx ()) + ii]);

      integral_h  += .25 * .25*area*dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
      integral_Ux += .25 * .25*area*dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
      integral_Uy += .25 * .25*area*dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
    }


    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){

        soll [ordh_nodal     (quadrant->gt (ii))] = std::max( soll [ordh_nodal     (quadrant->gt (ii))], inh );

        sol [ordh_nodal     (quadrant->gt (ii))] += integral_h;
        sol [ordUx_nodal    (quadrant->gt (ii))] += integral_Ux;
        sol [ordUy_nodal    (quadrant->gt (ii))] += integral_Uy;
      }

      else
      {
        // touch parent nodes to set up distributed vector structure
        sol [ordh_nodal   (quadrant->gparent(0,ii))] += integral_h*.5;
        sol [ordh_nodal   (quadrant->gparent(1,ii))] += integral_h*.5;

        sol [ordUx_nodal  (quadrant->gparent(0,ii))] += integral_Ux*.5;
        sol [ordUx_nodal  (quadrant->gparent(1,ii))] += integral_Ux*.5;

        sol [ordUy_nodal  (quadrant->gparent(0,ii))] += integral_Uy*.5;
        sol [ordUy_nodal  (quadrant->gparent(1,ii))] += integral_Uy*.5;
      }
    }
  }

  for (auto kk = 0; kk < mass.get_owned_data ().size (); kk++)
  {
    sol.get_owned_data ()[kk] /= mass.get_owned_data ()[kk];
  }

  // bim2a_solution_with_ghosts in quad_operators.cpp
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh_nodal,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx_nodal, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy_nodal);
  TOC("get nodal sol.");
  


  if (is_initial_refinement)
  {


    TIC();
    Q1 only_h (ln_nodes);
    bim2a_solution_with_ghosts (tmsh, only_h);
    for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
    {
      only_h(idx) = sol(ordh_nodal(idx));
    }
    only_h.assemble (replace_op);
    TOC("get separated sol.");


   
    TIC();
    auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);
    //q2_vec h_star = bim2c_quadtree_pde_recovered_solution (tmsh, only_h, dh);
    TOC ("gradient and hstar");


    TIC();
    // auto estimator = [& h_star, & only_h] (tmesh::quadrant_iterator q)
    // {
    //   return estimator_sol (q, h_star, only_h);
    // };
    auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
    {
      return estimator_grad(q, dh, only_h);
    };
    auto estimator_flux = [& only_h] (tmesh::quadrant_iterator q)
    {

      std::array<double,4> h_mesh = {0,0,0,0};
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! q->is_hanging (ii)){
          h_mesh[ii] = only_h[q->gt (ii)];
        }
        else
        {
          h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        }
      }

      const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

      return (basin_check); 
    };


    auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
    {

      std::array<double,4> h_mesh = {0,0,0,0};
      for (int ii = 0; ii < 4; ++ii)
      {
        if (! q->is_hanging (ii)){
          h_mesh[ii] = only_h[q->gt (ii)];
        }
        else
        {
          h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
        }
      }

      const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])==0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

      return (basin_check); 
    };


    tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, 1e-5, 4, 0, 0);
    //tmsh.set_metrics_marker (estimator, 1e-5, 4, 3, 1);
    tmsh.metrics_refine (1e7);  // RAFFINAMENTO (arg is max element)

    // tmsh.set_coarsen_marker (coarsen_function);
    // tmsh.set_refine_marker  (refine_function);
    // tmsh.coarsen (recursive, 1, 0);
    // tmsh.refine  (recursive, 1);
    TOC ("refine");
  
    // Ottengo i parametri della mesh corrente
    TIC();
    gn_nodes    = tmsh.num_global_nodes ();
    ln_nodes    = tmsh.num_owned_nodes ();
    ln_elements = tmsh.num_local_quadrants ();
    gn_elements = tmsh.num_global_quadrants ();
    TOC ("Obtaining new parameters");
  
  
    // Interpolo sol sulla nuova mesh
    TIC();
  
    Q0 sol_onehalf_ (ln_elements * 3);
    sol_onehalf_.assign (sol_onehalf_.size(), 0.0);

    distributed_vector dg_coefficients_ (ln_elements * 4 * 3), incr_dg_coefficients_(ln_elements * 4 * 3);
    dg_coefficients_.get_owned_data ().assign(dg_coefficients_.size(), 0.);
    incr_dg_coefficients_.get_owned_data ().assign(incr_dg_coefficients_.size(), 0.);


    Q1 sol_ (ln_nodes * 3);
    Q1 Z_ (ln_nodes);
    Q0 slope_x_(ln_elements);
    Q0 slope_y_(ln_elements);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {

      const auto iMe = elemental_mass_matrix_inverted(quadrant);

      const auto Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
      const auto Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
      const auto area = Dx * Dy;

      double xx_c=quadrant->centroid(0);
      double yy_c=quadrant->centroid(1); 

      slope_x_[quadrant->get_forest_quad_idx ()] = dem_slope_x[global_coord_2_raster(xx_c,yy_c)[0]];
      slope_y_[quadrant->get_forest_quad_idx ()] = dem_slope_y[global_coord_2_raster(xx_c,yy_c)[0]];

      std::array<double,4> dg_coeff_h = {0., 0., 0., 0.}, dg_coeff_Ux = {0., 0., 0., 0.}, dg_coeff_Uy = {0., 0., 0., 0.};
      for (int kk = 0; kk < 4; ++kk)
      {
        for (int ii = 0; ii < 4; ++ii)
        {
          double xx=quadrant->p(0,ii);
          double yy=quadrant->p(1,ii); 

          // trapezoidal rule
          if (! quadrant->is_hanging (ii)){
            dg_coeff_h [kk] += .25 * area * h0_fun   (xx, yy) * 1. * (ii == kk);
            dg_coeff_Ux[kk] += .25 * area * Ux0_fun  (xx, yy) * 1. * (ii == kk);
            dg_coeff_Uy[kk] += .25 * area * Uy0_fun  (xx, yy) * 1. * (ii == kk);
          }
          else
          {
            dg_coeff_h [kk] += .25 * area * h0_fun   (xx, yy) * .5 * (ii == kk);
            dg_coeff_Ux[kk] += .25 * area * Ux0_fun  (xx, yy) * .5 * (ii == kk);
            dg_coeff_Uy[kk] += .25 * area * Uy0_fun  (xx, yy) * .5 * (ii == kk);
          }

        }
      }

      std::array<double,4> dg_coeff_h_, dg_coeff_Ux_, dg_coeff_Uy_;
      for (int kk = 0; kk < 4; ++kk)
      {
        const auto & iMe_kk = iMe[kk];

        // matrice consistente
        // dg_coeff_h_ [kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_h .begin(), 0.); 
        // dg_coeff_Ux_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_Ux.begin(), 0.); 
        // dg_coeff_Uy_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), dg_coeff_Uy.begin(), 0.); 


        // matrice lumpata
        const double sum_iMe_kk = std::accumulate(iMe_kk.begin(), iMe_kk.end(), 0.);

        dg_coeff_h_ [kk] = sum_iMe_kk*dg_coeff_h [kk]; 
        dg_coeff_Ux_[kk] = sum_iMe_kk*dg_coeff_Ux[kk]; 
        dg_coeff_Uy_[kk] = sum_iMe_kk*dg_coeff_Uy[kk]; 

      }



      for (int kk = 0; kk < 4; ++kk)
      {
        dg_coefficients_[ordh  (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_h_ [kk];
        dg_coefficients_[ordUx (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_Ux_[kk];
        dg_coefficients_[ordUy (quadrant->get_global_quad_idx ()) + kk] = dg_coeff_Uy_[kk];

        incr_dg_coefficients_[ordh  (quadrant->get_global_quad_idx ()) + kk] = 0.;
        incr_dg_coefficients_[ordUx (quadrant->get_global_quad_idx ()) + kk] = 0.;
        incr_dg_coefficients_[ordUy (quadrant->get_global_quad_idx ()) + kk] = 0.;
      }



      for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
       quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
      {
        for (int kk = 0; kk < 4; ++kk)
        {
          dg_coefficients_[ordh  (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
          dg_coefficients_[ordUx (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
          dg_coefficients_[ordUy (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;

          incr_dg_coefficients_[ordh  (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
          incr_dg_coefficients_[ordUx (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
          incr_dg_coefficients_[ordUy (quadrant_nei->get_global_quad_idx ()) + kk] += 0.;
        }
      }




      for (int ii = 0; ii < 4; ++ii)
      {
        if (! quadrant->is_hanging (ii)){
          double xx=quadrant->p(0,ii);
          double yy=quadrant->p(1,ii);

          sol_ [ordh_nodal     (quadrant->gt (ii))] = h0_fun  (xx, yy);
          sol_ [ordUx_nodal    (quadrant->gt (ii))] = Ux0_fun (xx, yy);
          sol_ [ordUy_nodal    (quadrant->gt (ii))] = Uy0_fun (xx, yy);

          Z_[quadrant->gt (ii)] = dem[global_coord_2_raster(xx,yy)[0]]; 
        }
        
        else
        {
          sol_ [ordh_nodal   (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordh_nodal   (quadrant->gparent(1,ii))] += 0.;

          sol_ [ordUx_nodal  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUx_nodal  (quadrant->gparent(1,ii))] += 0.;

          sol_ [ordUy_nodal  (quadrant->gparent(0,ii))] += 0.;
          sol_ [ordUy_nodal  (quadrant->gparent(1,ii))] += 0.;

          Z_[quadrant->gparent(0,ii)] += 0.;
          Z_[quadrant->gparent(1,ii)] += 0.;
        }
      }
    }

    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordh_nodal,  false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUx_nodal, false);
    bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUy_nodal);

    dg_coefficients_.remap();
    dg_coefficients_.assemble(replace_op);

    Z_.remap();
    Z_.assemble(replace_op);

    incr_dg_coefficients_.remap();
    incr_dg_coefficients_.assemble(replace_op);



    sol                  = sol_;
    dg_coefficients      = dg_coefficients_;
    incr_dg_coefficients = incr_dg_coefficients_;
    sol_onehalf          = sol_onehalf_;
    Z                    = Z_;
    slope_x              = slope_x_;
    slope_y              = slope_y_;
  
    TOC ("compute initial condition");
  }

//./main $PWD inputs/dem_second_test.octbin.gz inputs/mask_in_vladi.octbin.gz

  
  //Q1 sol_dyn              = sol;
  Q1 Z_dyn                = Z;
  Q0 slope_x_dyn          = slope_x;
  Q0 slope_y_dyn          = slope_y;
  Q0 sol_onehalf_dyn      = sol_onehalf;

  distributed_vector dg_coefficients_dyn        = dg_coefficients;
  distributed_vector dg_coefficients_old_dyn    = dg_coefficients;
  distributed_vector dg_coefficients_oldold_dyn = dg_coefficients;
  distributed_vector incr_dg_coefficients_dyn   = incr_dg_coefficients;



  
  TG2_scheme stp(dg_coefficients_dyn, 
                 dg_coefficients_old_dyn, 
                 dg_coefficients_oldold_dyn, 
                 incr_dg_coefficients_dyn, 
                 sol_onehalf_dyn, 
                 ordh, ordUx, ordUy, ordh_nodal, ordUx_nodal, ordUy_nodal, 
                 Z_dyn, 
                 slope_x_dyn,
                 slope_y_dyn,
                 DELTAT, h_min, is_non_reflBC, is_bed_friction, is_stress_tensor, grav,
                 density, turbulence_coeff, surface_pressure, bed_friction_angle_rad, fluid_viscosity, yield_shear_stress);
  
  
  // Save initial conditions
  str = std::string(SAVE_DIR) + "/results/swel_h_%4.4d"; 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, soll, ordh_nodal);

  str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d"; 
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol, ordh_nodal);

  str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol, ordUx_nodal); 

  str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0);
  tmsh.octbin_export (filename, sol, ordUy_nodal);
  
  str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
  strcpy(arr, str.c_str());
  sprintf(filename, arr, 0); 
  tmsh.octbin_export (filename, Z_dyn);
  
  // MPI_Barrier (MPI_COMM_WORLD);
  // if (rank == 0) { print_timing_report (); }
  // MPI_Finalize (); 
  // return 0;
  


  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));

  
  // Time loop
  double time      = 0.0;
  double time_old  = 0.0;
  double time_oldd = 0.0;

  
  stp.set_dt (DELTAT);
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
  {
    stp.compute_dt(quadrant);
  }
  double max_dt = REDCDT * stp.dt;
  MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&max_dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);
  stp.set_dt(max_dt);
  time_old  -= stp.dt;
  time_oldd -= 2*stp.dt;
  stp.set_times(time, time_old, time_oldd);
  
  
  
  if(rank==0) {
    full_time_vector.push_back (0.0);
    save_time_vector.push_back (0.0);
  }
  int count = 0;
  
  double savecount = 0.0, space_adapt_count = 0.0;
  
  if (rank == 0)
  {
    std::cout << "start loop" << std::endl;
  }
  
  
  while (time < T)
  {
    
    
    // Reset increment, and limiter terms
    TIC();
    incr_dg_coefficients_dyn.get_owned_data ().assign (incr_dg_coefficients_dyn.get_owned_data ().size (), 0.0);
    incr_dg_coefficients_dyn.assemble (replace_op);
    TOC("Reset");
    TIC();
    
  
    // compute time step, 
    stp.set_dt (DELTAT);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      //stp.compute_dt(quadrant);
    }
    max_dt = REDCDT * stp.dt;
    stp.set_dt(max_dt); // deltat max
    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.dt), 1, MPI_DOUBLE, MPI_MIN, tmsh.comm);
    
    
    // // time adaptivity
    // if (is_time_adaptivity)
    // {
    //   stp.nu_htot = 0.;
    //   for (auto quadrant = tmsh.begin_quadrant_sweep ();
    //     quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    //   {
    //     stp.compute_dt_adaptive(quadrant);
    //   }
    //   MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&stp.nu_htot), 1, MPI_DOUBLE, MPI_SUM, tmsh.comm);
    //   const double rho_h = stp.nu_htot/((stp.time-stp.timed)*(stp.time-stp.timed));

    //   //std::cout << stp.nu_htot << " " << rank << std::endl;

    //   const auto candidate_dt = (5e-3/(rho_h+1e-7))*std::sqrt(1./(stp.time-stp.timed));
    //   stp.set_dt( std::min(std::isnan(candidate_dt) ? max_dt : candidate_dt, max_dt) );
    // }


    if (stp.dt == 0 && rank == 0)
    {
      std::cout << "dt has gone to zero, sorry, STOP!" << std::endl;
      exit( -1. );
    }


    // check save with given frequency
    //stp.set_dt((savecount+stp.dt)/SAVEDT>1 ? (stp.dt - std::fmod(savecount+stp.dt,SAVEDT) - SAVEDT*(std::floor(savecount+stp.dt/SAVEDT)-1))-SAVEDT : stp.dt);


    time_oldd = time_old;
    time_old = time;
    time += stp.dt; 
    savecount += stp.dt;
    space_adapt_count += stp.dt;
    MPI_Bcast (static_cast<void*> (&time),              1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Bcast (static_cast<void*> (&savecount),         1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Bcast (static_cast<void*> (&space_adapt_count), 1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Barrier (tmsh.comm); // tmsh.comm = MPI_COMM_WORLD
    
    // Print current time
    if(rank==0) 
    {
      std::cout << "TIME = " << time << ", dt = " << stp.dt << std::endl;
      full_time_vector.push_back (time);
    }
    
    
    
    
    // first step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.first_step(quadrant);
    }
    

    stp.set_times(time, time_old, time_oldd);
    // soldd_dyn = sold_dyn;
    // sold_dyn  = sol_dyn;
    

    // second step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.second_step(quadrant);
    } 
    incr_dg_coefficients_dyn.assemble ();
    TOC("Compute step"); 
    
    

    TIC();

    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      const auto iMe = elemental_mass_matrix_inverted(quadrant);

      std::array<double,4> result_h, result_Ux, result_Uy;
      for (int kk = 0; kk < 4; ++kk)
      {
        //std::cout << incr_dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + kk] << std::endl;// " " << incr_dg_coefficients_dyn[ordUx   (quadrant->get_global_quad_idx ()) + kk] << " " << incr_dg_coefficients_dyn[ordUy   (quadrant->get_global_quad_idx ()) + kk] << std::endl;

        result_h [kk] = incr_dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + kk];
        result_Ux[kk] = incr_dg_coefficients_dyn[ordUx  (quadrant->get_global_quad_idx ()) + kk];
        result_Uy[kk] = incr_dg_coefficients_dyn[ordUy  (quadrant->get_global_quad_idx ()) + kk];
      }

      std::array<double,4> result_h_, result_Ux_, result_Uy_;
      for (int kk = 0; kk < 4; ++kk)
      {
        const auto & iMe_kk = iMe[kk];

        result_h_ [kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), result_h .begin(), 0.); 
        result_Ux_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), result_Ux.begin(), 0.); 
        result_Uy_[kk] = std::inner_product(iMe_kk.begin(), iMe_kk.end(), result_Uy.begin(), 0.); 

        // // matrice lumpata
        // const double sum_iMe_kk = std::accumulate(iMe_kk.begin(), iMe_kk.end(), 0.);
        
        // result_h_ [kk] = sum_iMe_kk*result_h [kk]; 
        // result_Ux_[kk] = sum_iMe_kk*result_Ux[kk]; 
        // result_Uy_[kk] = sum_iMe_kk*result_Uy[kk]; 
      }

      for (int kk = 0; kk < 4; ++kk)
      {
        dg_coefficients_dyn[ordh  (quadrant->get_global_quad_idx ()) + kk] += stp.dt * result_h_ [kk];
        dg_coefficients_dyn[ordUx (quadrant->get_global_quad_idx ()) + kk] += stp.dt * result_Ux_[kk];
        dg_coefficients_dyn[ordUy (quadrant->get_global_quad_idx ()) + kk] += stp.dt * result_Uy_[kk];
      }

    }

    
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
      {
        //std::cout << dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + ii] << " " << dg_coefficients_dyn[ordUx   (quadrant->get_global_quad_idx ()) + ii] << " " << dg_coefficients_dyn[ordUy   (quadrant->get_global_quad_idx ()) + ii] << std::endl;

        //if (!std::isnan(dg_coefficients_dyn [ordh (quadrant->get_global_quad_idx ()) + ii]))
        //std::cout << "here!! " << dg_coefficients_dyn [ordh (quadrant->get_global_quad_idx ()) + ii] << std::endl;
        dg_coefficients_dyn [ordh (quadrant->get_global_quad_idx ()) + ii] *= (dg_coefficients_dyn [ordh (quadrant->get_global_quad_idx ()) + ii]>0);
      }
    }
    dg_coefficients_dyn.assemble (replace_op);
    TOC("Apply increment");

    //exit(1);
    
    
    // Save solution
    if (savecount >= SAVEDT || time>=0.5) 
    {
      TIC();
      if (rank == 0)
        std::cout << "savecount = " << savecount << std::endl;
      count++;
      save_time_vector.push_back (time);

      TIC();
      Q1 sol_(ln_nodes * 3);
      sol_.get_owned_data ().assign(sol_.size(), 0.);

      Q1 soll_(ln_nodes * 3);
      soll_.get_owned_data ().assign(soll_.size(), 0.);

      Q1 mass (ln_nodes * 3);
      bim2a_mass_vector (tmsh, mass, ordh_nodal );
      bim2a_mass_vector (tmsh, mass, ordUx_nodal);
      bim2a_mass_vector (tmsh, mass, ordUy_nodal);
      mass.assemble ();


      for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
      {

        std::array<double,4> isdof_or_hanging;

        const auto Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
        const auto Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
        const auto area = Dx * Dy;

        for (int ii = 0; ii < 4; ++ii){
          if (! quadrant->is_hanging (ii)){
            isdof_or_hanging[ii] = 1.;
          } else {
            isdof_or_hanging[ii] = .5;
          }
        }

        // midpoint integration
        double inh = -1e4, inUx = -1e4, inUy = -1e4;
        double integral_h = 0., integral_Ux = 0., integral_Uy = 0.;
        for (int ii = 0; ii < 4; ++ii)
        {
          //if (dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + ii] >0)
          //std::cout << ii << " " << integral_h << " " << dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + ii] << std::endl;

          inh  = std::max(inh,  dg_coefficients_dyn[ordh    (quadrant->get_global_quad_idx ()) + ii]);
          inUx = std::max(inUx, dg_coefficients_dyn[ordUx   (quadrant->get_global_quad_idx ()) + ii]);
          inUy = std::max(inUy, dg_coefficients_dyn[ordUy   (quadrant->get_global_quad_idx ()) + ii]);

          integral_h  += .25 * .25*area*dg_coefficients_dyn[ordh   (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
          integral_Ux += .25 * .25*area*dg_coefficients_dyn[ordUx  (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
          integral_Uy += .25 * .25*area*dg_coefficients_dyn[ordUy  (quadrant->get_global_quad_idx ()) + ii]*isdof_or_hanging[ii];
        }

        

        for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){

            soll_ [ordh_nodal     (quadrant->gt (ii))] = std::max( soll_ [ordh_nodal     (quadrant->gt (ii))], inh );

            sol_ [ordh_nodal     (quadrant->gt (ii))] += integral_h;
            sol_ [ordUx_nodal    (quadrant->gt (ii))] += integral_Ux;
            sol_ [ordUy_nodal    (quadrant->gt (ii))] += integral_Uy;
          }

          else
          {
             
            // touch parent nodes to set up distributed vector structure
            sol_ [ordh_nodal   (quadrant->gparent(0,ii))] += integral_h*.5;
            sol_ [ordh_nodal   (quadrant->gparent(1,ii))] += integral_h*.5;

            sol_ [ordUx_nodal  (quadrant->gparent(0,ii))] += integral_Ux*.5;
            sol_ [ordUx_nodal  (quadrant->gparent(1,ii))] += integral_Ux*.5;

            sol_ [ordUy_nodal  (quadrant->gparent(0,ii))] += integral_Uy*.5;
            sol_ [ordUy_nodal  (quadrant->gparent(1,ii))] += integral_Uy*.5;
          }
        }
      }

      for (auto kk = 0; kk < mass.get_owned_data ().size (); kk++)
      {
        sol_.get_owned_data ()[kk] /= mass.get_owned_data ()[kk];
      }

      // bim2a_solution_with_ghosts in quad_operators.cpp
      bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordh_nodal,  false);
      bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUx_nodal, false);
      bim2a_solution_with_ghosts (tmsh, sol_, replace_op, ordUy_nodal);
      TOC("get nodal sol.");


      str = std::string(SAVE_DIR) + "/results/swel_h_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,   count);
      tmsh.octbin_export (filename, soll_, ordh_nodal);


     
      str = std::string(SAVE_DIR) + "/results/swe_h_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,   count);
      tmsh.octbin_export (filename, sol_, ordh_nodal);
      
      str = std::string(SAVE_DIR) + "/results/swe_Ux_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_, ordUx_nodal); 
      
      str = std::string(SAVE_DIR) + "/results/swe_Uy_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, sol_, ordUy_nodal);
      
      str = std::string(SAVE_DIR) + "/results/swe_Z_%4.4d";
      strcpy(arr, str.c_str());
      sprintf(filename, arr,  count);
      tmsh.octbin_export (filename, Z_dyn);

      savecount = 0.0;
      TOC("Exporting solution");
      
      //exit(1);

    }

    /*
    if (is_space_adaptivity && space_adapt_count >= SPACE_ADAPTDT)
    {

      TIC();
      Q1 only_h (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_h);
      for (auto idx = only_h.get_range_start (); idx != only_h.get_range_end (); ++idx)
      {
        only_h(idx) = sol_dyn(ordh(idx));
      }
      only_h.assemble (replace_op);

      Q1 only_Ux (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Ux);
      for (auto idx = only_Ux.get_range_start (); idx != only_Ux.get_range_end (); ++idx)
      {
        only_Ux(idx) = sol_dyn(ordUx(idx));
      }
      only_Ux.assemble (replace_op);

      Q1 only_Uy (ln_nodes);
      bim2a_solution_with_ghosts (tmsh, only_Uy);
      for (auto idx = only_Uy.get_range_start (); idx != only_Uy.get_range_end (); ++idx)
      {
        only_Uy(idx) = sol_dyn(ordUy(idx));
      }
      only_Uy.assemble (replace_op);      
      TOC("get separated sol.");

      TIC();
      std::set<int> global_index_quad;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep ();
           ++q)
      {
        quadrant_marker_list(q, only_h,  only_Ux, only_Uy, stp.dt, global_index_quad);
      }      
      TOC("front track.");  
      
      
      TIC();
      auto dh = bim2c_quadtree_pde_recovered_gradient (tmsh, only_h);
      //q2_vec h_star = bim2c_quadtree_pde_recovered_solution (tmsh, only_h, dh);
      TOC ("gradient and hstar");
      
      TIC();
      // auto estimator = [& h_star, & only_h] (tmesh::quadrant_iterator q)
      // {
      //   return estimator_sol (q, h_star, only_h);
      // };
      auto estimator = [& dh, & only_h] (tmesh::quadrant_iterator q)
      {
        return estimator_grad(q, dh, only_h);
      };

      // auto estimator_flux = [& h_star, & only_h] (tmesh::quadrant_iterator q)
      // {
      //   std::array<double,4> h_mesh = {0,0,0,0};
      //   for (int ii = 0; ii < 4; ++ii)
      //   {
      //     if (! q->is_hanging (ii)){
      //       h_mesh[ii] = only_h[q->gt (ii)];
      //     }
      //     else
      //     {
      //       h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
      //     } 
      //   }

      //   const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 
      //   //const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])>h_min && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])<=h_min) ? 1 : 0; 

      //   return (basin_check); 
      // };


      auto estimator_flux = [& global_index_quad] (tmesh::quadrant_iterator q)
      {
        if ( global_index_quad.find(q->get_global_quad_idx ()) != global_index_quad.end() )
        {
          return 1;
        }
        return 0;
      };


      auto dry_function = [& only_h] (tmesh::quadrant_iterator q)
      {

        std::array<double,4> h_mesh = {0,0,0,0};
        for (int ii = 0; ii < 4; ++ii)
        {
          if (! q->is_hanging (ii)){
            h_mesh[ii] = only_h[q->gt (ii)];
          }
          else
          {
            h_mesh[ii] = .5 * ( only_h[q->gparent(0,ii)] + only_h[q->gparent(1,ii)] );
          }
        }

        const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])<h_min && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])<h_min) ? 1 : 0; 
        //const auto basin_check = ((h_mesh[0]+h_mesh[1]+h_mesh[2]+h_mesh[3])==0 && (h_mesh[0]*h_mesh[1]*h_mesh[2]*h_mesh[3])==0) ? 1 : 0; 

        return (basin_check); 
      };


      tmsh.set_metrics_marker_flux_lim (estimator, estimator_flux, dry_function, mesh_size_dry, mesh_size_wet, mesh_size_interface, 1e-5, 4, 0, 0);
      //tmsh.set_metrics_marker (estimator, 1e-5, 4, 3, 1); 
      tmsh.metrics_refine (1e7);  // RAFFINAMENTO (arg is max element)

      // tmsh.set_coarsen_marker (coarsen_function);
      // tmsh.set_refine_marker  (refine_function);
      // tmsh.coarsen (recursive, 1, 0);
      // tmsh.refine  (recursive, 1);
      TOC ("refine");

      // Ottengo i parametri della mesh corrente
      TIC();
      gn_nodes    = tmsh.num_global_nodes ();
      ln_nodes    = tmsh.num_owned_nodes ();
      ln_elements = tmsh.num_local_quadrants ();
      gn_elements = tmsh.num_global_quadrants ();
      TOC ("Obtaining new parameters");
      
      
      // Interpolo sol sulla nuova mesh
      TIC();
      Q1 sol (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
      interpolate_vector (tmsh, sol_dyn, sol, ordh);
      interpolate_vector (tmsh, sol_dyn, sol, ordUx);
      interpolate_vector (tmsh, sol_dyn, sol, ordUy);
      //sol.assemble (replace_op);
      
      Q1 sold (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, sold, replace_op, ordUy);
      interpolate_vector (tmsh, sold_dyn, sold, ordh);
      interpolate_vector (tmsh, sold_dyn, sold, ordUx);
      interpolate_vector (tmsh, sold_dyn, sold, ordUy);
      //sold.assemble (replace_op);
      
      
      Q1 soldd (ln_nodes * 3);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, soldd, replace_op, ordUy);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordh);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUy);
      interpolate_vector (tmsh, soldd_dyn, soldd, ordUx);
      //soldd.assemble (replace_op);
      
      
      Q1 incr (ln_nodes * 3);
      incr.get_owned_data ().assign (incr.get_owned_data ().size(), 0.0);
      incr.assemble ();

      std::vector<std::array<double,4>> incr_anti_diff (ln_elements * 3);

      Q1 P_plus (ln_nodes * 3);
      P_plus.get_owned_data ().assign (P_plus.get_owned_data ().size(), 0.0);
      P_plus.assemble ();  
      
      Q1 P_minus (ln_nodes * 3);
      P_minus.get_owned_data ().assign (P_minus.get_owned_data ().size(), 0.0);
      P_minus.assemble ();      
      
      Q1 mass (ln_nodes * 3);
      bim2a_mass_vector (tmsh, mass, ordh);
      bim2a_mass_vector (tmsh, mass, ordUx);
      bim2a_mass_vector (tmsh, mass, ordUy);
      mass.assemble ();
      
      Q0 sol_onehalf (ln_elements * 3);
      sol_onehalf.assign (sol_onehalf.size(), 0.0);


      Q1 Z (ln_nodes);
      Q0 slope_x(ln_elements);
      Q0 slope_y(ln_elements);
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
      {

        double xx_c=quadrant->centroid(0);
        double yy_c=quadrant->centroid(1); 

        slope_x[quadrant->get_forest_quad_idx ()] = dem_slope_x[global_coord_2_raster(xx_c,yy_c)[0]];
        slope_y[quadrant->get_forest_quad_idx ()] = dem_slope_y[global_coord_2_raster(xx_c,yy_c)[0]];
        
        for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            double xx=quadrant->p(0,ii);
            double yy=quadrant->p(1,ii);
            Z[quadrant->gt (ii)] = dem[global_coord_2_raster(xx,yy)[0]]; 
          }
           
          else
          {
            Z[quadrant->gparent(0,ii)] += 0.;
            Z[quadrant->gparent(1,ii)] += 0.;
          }
        }
      }
      
      
      
      bim2a_solution_with_ghosts (tmsh, Z, replace_op);
      
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);

      bim2a_solution_with_ghosts (tmsh, P_plus, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, P_plus, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, P_plus, replace_op, ordUy);

      bim2a_solution_with_ghosts (tmsh, P_minus, replace_op, ordh,  false);
      bim2a_solution_with_ghosts (tmsh, P_minus, replace_op, ordUx, false);
      bim2a_solution_with_ghosts (tmsh, P_minus, replace_op, ordUy);
      
       
      sol_dyn              = sol;
      sold_dyn             = sold;
      soldd_dyn            = soldd;
      incr_dyn             = incr;
      incr_anti_diff_dyn   = incr_anti_diff;
      P_plus_dyn           = P_plus;
      P_minus_dyn          = P_minus;
      mass_dyn             = mass;
      sol_onehalf_dyn      = sol_onehalf;
      Z_dyn                = Z;
      slope_x_dyn          = slope_x;
      slope_y_dyn          = slope_y;

      space_adapt_count = 0.0;
      
      TOC ("Interpolation");
      
    }*/
      
    
    
    
  }
  
  
  if (rank == 0)
  {
    str = std::string(SAVE_DIR) + "/results/timesteps.octbin"; 
    strcpy(arr, str.c_str());
    sprintf(filename, arr, 0);

    ColumnVector vtmp (save_time_vector.size ());
    std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
    octave_io_mode m;
    octave_io_open (filename, gz_write_mode, &m);
    octave_save ("save_time", vtmp);
    vtmp.resize (full_time_vector.size ());
    std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
    octave_save ("full_time", vtmp);
    octave_io_close ();
  }
  
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  MPI_Finalize ();
  return 0;
  
}






