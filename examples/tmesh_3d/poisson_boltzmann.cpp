#include <mpi.h>

#include "pb_class.h"
#include "pqr_parser.cpp"

static char filename[255];

double crossings_t::start = 0.;
double crossings_t::end = 0.;

std::vector<NS::Atom> a;
NS::NanoShaper crossings_t::ns(a, NS::skin, 0.45, 2., 1);

int_coord_t ray_cache_t::count_cache = 0;
int_coord_t ray_cache_t::count_new = 0;

void
crossings_to_json (nlohmann::json& j, const crossings_t& cr_t);
void 
map_to_json(nlohmann::json& j, const std::map<std::array<double, 2> , crossings_t, map_compare>& r);

void
print_map(const std::map<std::array<double, 2> , crossings_t, map_compare>& r); 
   
int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
 
  poisson_boltzmann pb;
  ray_cache_t ray_cache;
  
  if (pb.parse_options (argc, argv))
    return 1;
  
  std::ifstream inputfile (pb.pqrfilename);
  read_atoms_from_pqr (inputfile, pb.atoms);
  inputfile.close ();
  
  if (rank == 0)
  {
    std::cout << "Atom : " << std::endl;
    write_atoms_to_pqr (std::cout, pb.atoms);
    pb.print_options ();
  }
  
  MPI_Barrier (mpicomm);
    
  TIC ();
  if (pb.mesh_shape == 1)
     pb.create_mesh ();
  else if (pb.mesh_shape == 0)
     pb.create_cubic_mesh ();
  else 
  {
     std::cerr << "Invalid mesh shape selected" << std::endl;
     return 1;
  }
  TOC ("create_mesh");

  TIC ();
  pb.init_tmesh ();
  TOC ("init_tmesh");
  
  crossings_t::start = pb.l_c[1]/2.;
  crossings_t::end = pb.r_c[1]/2.;
 
  NS::NanoShaper ns2 (pb.atoms, pb.surf_type, pb.skin_param, pb.stern_layer, pb.numberOfThreads);
  crossings_t::ns = ns2;
  crossings_t::ns.setConfig<double>("Grid_scale", 3.0 );
  crossings_t::ns.buildAnalyticalSurface();
  unsigned y_direction = 1; 
  crossings_t::ns.setDirection(y_direction);
  std::cout << "\n" << std::endl;
  
  MPI_Barrier (mpicomm);
  
  TIC ();
  //pb.refine_surface ();
  pb.refine_surface_ns (ray_cache);
  TOC ("refine surface");

  TIC ();
  //pb.create_markers ();
  pb.create_markers_ns (ray_cache);
  TOC ("create element markers");
  
  TIC ();
  //pb.export_ls_tmesh ();
  pb.export_ns_tmesh (ray_cache);
  TOC ("export ls tmesh");

  TIC ();
  pb.export_marked_tmesh ();
  TOC ("export marked tmesh");
  
  TIC ();
  if (pb.linear_solver_name == "mumps")
     pb.mumps_compute_electric_potential ();
  else if (pb.linear_solver_name == "lis")
     pb.lis_compute_electric_potential ();
  else 
  {
     std::cerr << "Invalid linear solver selected" << std::endl;
     return 1;
  }
  TOC ("compute electric potential");
  
  TIC ();
  pb.export_p4est ();
  TOC ("export p4est");
    
  if (rank == 0) { print_timing_report(); }
  MPI_Barrier (mpicomm);
     
  MPI_Finalize ();
  
  /*
  //Save ray_cache:
  nlohmann::json j;
  map_to_json (j, ray_cache.rays);
  //std::cout << "Cache size: " << pb.ray_cache.rays.size() << std::endl;
  //std::cout << j << std::endl;
  
  std::ofstream ray_cached_file;
  ray_cached_file.open ("ray_cache.txt");
  
  if (ray_cached_file.is_open ())
    ray_cached_file << j;
    
  ray_cached_file.close ();*/
  
  print_map (ray_cache.rays); 
  
  return 0;
  
}

void
print_map(const std::map<std::array<double, 2> , crossings_t, map_compare>& r) 
{
  std::ofstream ray_cached_file;
  ray_cached_file.open ("ray_cache_2.txt");
  if (ray_cached_file.is_open ())
  {
    ray_cached_file << "Count cached rays: " << ray_cache_t::count_cache << std::endl;
    ray_cached_file << "Count new rays: " << ray_cache_t::count_new << std::endl;
    for (auto it : r)
    {
      ray_cached_file << "[[" << it.first.at(0) << ", " << it.first.at(1) << "]";
      for (int i = 0; i < it.second.inters.size (); i++) 
        ray_cached_file << ", " << it.second.inters[i].first;
      ray_cached_file << "]" << std::endl;
    }
  }
  ray_cached_file.close (); 
}

void
crossings_to_json (nlohmann::json& j, const crossings_t& cr_t) 
{
  if (cr_t.inters.begin () != cr_t.inters.end ()) //nota: non sta mettendo tutto il vettore inters, ma solo il primo valore. non va a capo. 
    j += nlohmann::json{{"flags", cr_t.flags}, {"intersections", cr_t.inters.begin()->first}};
  else
    j += nlohmann::json{{"flags", cr_t.flags}, {"intersections", cr_t.flags}};
}

void 
map_to_json(nlohmann::json& j, const std::map<std::array<double, 2> , crossings_t, map_compare>& r) 
{
  for (auto it : r)
  {
  j += nlohmann::json{{it.first}};
  crossings_to_json(j, it.second);
  }
  
  /*int count = 0;
  for (auto it : r)
  {
    count ++;
    auto i = it.second.inters.cbegin();
    
    j += nlohmann::json{{"ray", it.first}, {"inters", (i->first)}};
  }*/
}
