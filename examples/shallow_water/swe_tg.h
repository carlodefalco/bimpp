// Setting parameters
constexpr double pi = 3.14159265358979323846264338327950288;
constexpr int NUM_REFINEMENTS = 4;
constexpr double SAVEDT  = 10;
constexpr double DELTAT =  1.;
constexpr double REDCDT =  .5;
constexpr double T      =  1000;


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

const dirichlet_bcs bcsh;
const dirichlet_bcs bcsUx = {{0, 0, [] (double, double) {return 0;}},
                             {9, 1, [] (double, double) {return 0;}}};

const dirichlet_bcs bcsUy= {{0, 3, [] (double, double) {return 0;}},
                            {0, 2, [] (double, double) {return 0;}},
                            {1, 3, [] (double, double) {return 0;}},
                            {1, 2, [] (double, double) {return 0;}},
                            {2, 3, [] (double, double) {return 0;}},
                            {2, 2, [] (double, double) {return 0;}},
                            {3, 3, [] (double, double) {return 0;}},
                            {3, 2, [] (double, double) {return 0;}},
                            {4, 3, [] (double, double) {return 0;}},
                            {4, 2, [] (double, double) {return 0;}},
                            {5, 3, [] (double, double) {return 0;}},
                            {5, 2, [] (double, double) {return 0;}},
                            {6, 3, [] (double, double) {return 0;}},
                            {6, 2, [] (double, double) {return 0;}},
                            {7, 3, [] (double, double) {return 0;}},
                            {7, 2, [] (double, double) {return 0;}},
                            {8, 3, [] (double, double) {return 0;}},
                            {8, 2, [] (double, double) {return 0;}},
                            {9, 3, [] (double, double) {return 0;}},
                            {9, 2, [] (double, double) {return 0;}}};

//double Z_fun  (double xx, double yy)  { return (3. - 3. * std::sin (pi * xx / 2. / 100.)); }
double Z_fun  (double xx, double yy)  { return (0.); }
double h0_fun (double xx, double yy)  { return std::max (0., (8. - std::sin (pi * xx / 2. / 400.) - Z_fun (xx, yy))); }
double Ux0_fun (double xx, double yy) { return 0.; }
double Uy0_fun (double xx, double yy) { return 0.; }

