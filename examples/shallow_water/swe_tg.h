// Setting parameters
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr int NUM_REFINEMENTS = 4;
constexpr double KAPPA  = -0.01;
constexpr double DELTAT =  0.005;
constexpr double T      = 2;


// Connectivity of local element
constexpr p4est_topidx_t simple_conn_num_vertices = 22;
constexpr p4est_topidx_t simple_conn_num_trees = 10;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {-400,  0,
   -400, 80,
   -320,  0,
   -320, 80,
   -240,  0,
   -240, 80,
   -160,  0,
   -160, 80,
    -80,  0,
    -80, 80,
      0,  0,
      0, 80,
     80,  0,
     80, 80,
    160,  0,
    160, 80,
    240,  0,
    240, 80,
    320,  0,
    320, 80,
    400,  0,
    400, 80};

const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {  1,    3,    4,    2,    1,
     3,    5,    6,    4,    1,
     5,    7,    8,    6,    1,
     7,    9,   10,    8,    1,
     9,   11,   12,   10,    1,
     11,   13,   14,   12,    1,
     13,   15,   16,   14,    1,
     15,   17,   18,   16,    1,
     17,   19,   20,   18,    1,
     19,   21,   22,   20,    1};


using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>;         // Typedef for local q_0 vector

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

/*
static
std::array<double, 4> rhsloc;

void
loop_over_vector (tmesh& mesh,
                  const Q0& f,
                  const Q1& g,
                  Q1& rhs,
                  const ordering& ord)
{
  rhsloc.fill (0.0);

  double f_loc = 0;
  std::array<double, 4> g_loc;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      f_loc = f[quadrant->get_forest_quad_idx ()];

      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            g_loc[n] = g[quadrant->gt (n)];
          else
            g_loc[n] = 0.5 * (g[quadrant->gparent (0, n)] +
                              g[quadrant->gparent (1, n)]);
        }

      //bim2a_rhs_loc (quadrant, f_loc, g_loc, rhsloc);
      assemble_vector (quadrant, rhsloc, rhs, ord);
    }
}
*/
