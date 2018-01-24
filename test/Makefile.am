BIMMUMPSLIB = ../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

CLEANFILES = 


check_PROGRAMS = io_test bgs_test sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 bim3_stationary_test \
                 bim3_stationary_test2 bim3_stationary_test2_osc bim3_stationary_test3 \
                 bim3_stationary_test3_osc bim3_stationary_test4 bim3_timedependent_test \
                 bim3_timedependent_test2 bim3_timedependent_test2_osc bim3_timedependent_test3 \
                 bim3_timedependent_test3_osc 

if HAVE_LIS
  check_PROGRAMS += lis_stationary_test lis_timedependent_test linear_solver_test \
                    linear_solver_test_repeated nonlinear_solver_test nonlinear_solver_test2 \
                    nonlinear_solver_test3 nonlinear_solver_test4 nonlinear_solver_eq_test
endif

if HAVE_P4EST

  check_PROGRAMS += p4est_quad_index_test p4est_conn_test_hanging p4est_hanging_test_2 p4est_operator_test \
                    p4est_conn_test p4est_hanging_test p4est_refine_segment_test p4est_coarsen_test \
                    p4est_iterator_test p4est_ref_test  p4est_load_save_test p4est_edge_id_test \
                    p4est_operator_test_2 p4est_operator_test_3 p4est_operator_segment_test \
                    p4est_neighbor_iterator_test p4est_estimator_test_1 p4est_up_test \
                    p4est_ref_lnodes_test p4est_test_discontinuous p4est_estimator_test_1_uniform \
                    p4est_metrics_test_1 \
                    p4est_dr_test_1_marker p4est_dr_test_1_metrics \
                    p4est_dr_test_2_marker p4est_dr_test_2_metrics \
                    p4est_dr_test_3_marker p4est_dr_test_3_metrics \
                    p4est_dr_test_4_marker p4est_dr_test_4_metrics \
                    p4est_adr_test_1_marker p4est_adr_test_1_metrics \
                    p4est_adr_test_2_marker p4est_adr_test_2_metrics

endif

# schur_test_SOURCES = schur_test.cpp
# CLEANFILES += 
# schur_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB)

io_test_SOURCES =  io_test.cpp
io_test_CPPFLAGS =  $(AM_CPPFLAGS)
io_test_LDADD = $(BIMLIB) $(BIMIOLIB)
CLEANFILES += tmp_io_test.octbin.gz

bgs_test_SOURCES = bgs_test.cpp
CLEANFILES += 
bgs_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) 

sparse_test_SOURCES = sparse_test.cpp
CLEANFILES += first_solve.txt second_solve.txt
sparse_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG.m mesh_out.m

bim3_test2_SOURCES = bim3_test2.cpp
bim3_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG2.m mesh_out.m

bim3_test3_SOURCES = bim3_test3.cpp
bim3_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG3.m mesh_out.m

bim3_test4_SOURCES = bim3_test4.cpp
bim3_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += UP.m mesh_out.m

bim3_stationary_test_SOURCES = bim3_stationary_test.cpp
bim3_stationary_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution.txt

bim3_stationary_test2_SOURCES = bim3_stationary_test2.cpp
bim3_stationary_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution2.txt mesh_out_cube.m

bim3_stationary_test2_osc_SOURCES = bim3_stationary_test2_osc.cpp
bim3_stationary_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution2_OSC.txt

bim3_stationary_test3_SOURCES = bim3_stationary_test3.cpp
bim3_stationary_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution3.txt

bim3_stationary_test3_osc_SOURCES = bim3_stationary_test3_osc.cpp
bim3_stationary_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution3_OSC.txt

bim3_stationary_test4_SOURCES = bim3_stationary_test4.cpp
bim3_stationary_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution4.txt

bim3_timedependent_test_SOURCES = bim3_timedependent_test.cpp
bim3_timedependent_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest1.txt

bim3_timedependent_test2_SOURCES = bim3_timedependent_test2.cpp
bim3_timedependent_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest2.txt

bim3_timedependent_test2_osc_SOURCES = bim3_timedependent_test2_osc.cpp
bim3_timedependent_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest2_osc.txt

bim3_timedependent_test3_SOURCES = bim3_timedependent_test3.cpp
bim3_timedependent_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest3.txt

bim3_timedependent_test3_osc_SOURCES = bim3_timedependent_test3_osc.cpp
bim3_timedependent_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest3_osc.txt

if HAVE_LIS

  lis_stationary_test_SOURCES = lis_stationary_test.cpp
  lis_stationary_test_LDADD = $(BIMLIB) $(BIMLISLIB)
  CLEANFILES += LIS_solution.txt

  lis_timedependent_test_SOURCES = lis_timedependent_test.cpp
  lis_timedependent_test_LDADD = $(BIMLIB) $(BIMLISLIB)
  CLEANFILES += Lis_Solution_TimeTest.txt

  linear_solver_test_SOURCES = linear_solver_test.cpp
  linear_solver_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += MUMPS_solution.txt LIS_solution.txt

  linear_solver_test_repeated_SOURCES = linear_solver_test_repeated.cpp
  linear_solver_test_repeated_LDADD =  $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)

  nonlinear_solver_test_SOURCES =  nonlinear_solver_test.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test_LDADD = $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  CLEANFILES += NonLinear_Plaplacian_Test1.txt

  nonlinear_solver_test2_SOURCES =  nonlinear_solver_test2.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test2_LDADD = $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  CLEANFILES += NonLinear_Plaplacian_Test2.txt

  nonlinear_solver_test3_SOURCES =  nonlinear_solver_test3.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test3_LDADD = $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  CLEANFILES += NonLinear_Plaplacian_Test3.txt

  nonlinear_solver_test4_SOURCES =  nonlinear_solver_test4.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test4_LDADD = $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  CLEANFILES += NonLinear_Plaplacian_Test4.txt

  nonlinear_solver_eq_test_SOURCES =  nonlinear_solver_eq_test.cpp equation_class.h equation_class.cpp
  nonlinear_solver_eq_test_LDADD =  $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  CLEANFILES += NonLinear_Equation_Test.txt

endif

if HAVE_P4EST

  p4est_operator_segment_test_SOURCES =  p4est_operator_segment_test.cpp
  p4est_operator_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_operator_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_operator_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_operator_segment_test*.octbin.gz p4est_operator_segment_test.pvtu p4est_operator_segment_test.visit p4est_operator_segment_test*.vtu

  p4est_quad_index_test_SOURCES =  p4est_quad_index_test.cpp 
  p4est_quad_index_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_quad_index_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_quad_index_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_quad_index_test*.octbin.gz p4est_quad_index_test.pvtu p4est_quad_index_test.visit p4est_quad_index_test_*.vtu

  p4est_operator_test_3_SOURCES =  p4est_operator_test_3.cpp 
  p4est_operator_test_3_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_operator_test_3_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_operator_test_3_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_operator_test_3*.octbin.gz p4est_operator_test_3.pvtu p4est_operator_test_3.visit p4est_operator_test_3_*.vtu

  p4est_operator_test_2_SOURCES =  p4est_operator_test_2.cpp
  p4est_operator_test_2_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_operator_test_2_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_operator_test_2_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_operator_test_2*.octbin.gz p4est_operator_test_2.pvtu p4est_operator_test_2.visit p4est_operator_test_2_*.vtu

  p4est_conn_test_hanging_SOURCES =  p4est_conn_test_hanging.cpp 
  p4est_conn_test_hanging_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_conn_test_hanging_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_conn_test_hanging_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_conn_test_hanging*.octbin.gz p4est_conn_test_hanging.pvtu p4est_conn_test_hanging.visit p4est_conn_test_hanging_*.vtu

  p4est_hanging_test_2_SOURCES =  p4est_hanging_test_2.cpp
  p4est_hanging_test_2_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_hanging_test_2_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_hanging_test_2_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_hanging_test_2.octbin.gz p4est_hanging_test_2.pvtu p4est_hanging_test_2.visit p4est_hanging_test_2_*.vtu

  p4est_operator_test_SOURCES =  p4est_operator_test.cpp
  p4est_operator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_operator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_operator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_operator_test*.octbin.gz p4est_operator_test.pvtu p4est_operator_test.visit p4est_operator_test_*.vtu

  p4est_conn_test_SOURCES =  p4est_conn_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_conn_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_conn_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_conn_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_conn_test*.octbin.gz p4est_conn_test.pvtu p4est_conn_test.visit p4est_conn_test_*.vtu

  p4est_hanging_test_SOURCES =  p4est_hanging_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_hanging_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_hanging_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_hanging_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_hanging_test.octbin.gz p4est_hanging_test.pvtu p4est_hanging_test.visit p4est_hanging_test_*.vtu

  p4est_refine_segment_test_SOURCES =  p4est_refine_segment_test.cpp
  p4est_refine_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_refine_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_refine_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_unitsquare.octbin.gz p4est_refine_segment_test*.pvtu p4est_refine_segment_test*.visit p4est_refine_segment_test_*.vtu

  p4est_coarsen_test_SOURCES =  p4est_coarsen_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_coarsen_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_coarsen_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_coarsen_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_coarsen_test.octbin.gz p4est_coarsen_test.pvtu p4est_coarsen_test.visit p4est_coarsen_test_*.vtu

  p4est_iterator_test_SOURCES =  p4est_iterator_test.cpp 
  p4est_iterator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_iterator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_iterator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_iterator_test.octbin.gz p4est_iterator_test.pvtu p4est_iterator_test.visit p4est_iterator_test_*.vtu

  p4est_edge_id_test_SOURCES =  p4est_edge_id_test.cpp 
  p4est_edge_id_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_edge_id_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_edge_id_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_edge_id_test.octbin.gz p4est_edge_id_test.pvtu p4est_edge_id_test.visit p4est_edge_id_test_*.vtu

  p4est_ref_test_SOURCES =  p4est_ref_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_ref_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_ref_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_ref_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_ref_test.octbin.gz p4est_ref_test.pvtu p4est_ref_test.visit p4est_ref_test_*.vtu

  p4est_ref_lnodes_test_SOURCES =  p4est_ref_lnodes_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_ref_lnodes_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_ref_lnodes_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_ref_lnodes_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_ref_lnodes_test.octbin*.gz p4est_ref_lnodes_test*.pvtu p4est_ref_lnodes_test*.visit p4est_ref_lnodes_test_*.vtu

  p4est_load_save_test_SOURCES =  p4est_load_save_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_load_save_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_load_save_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_load_save_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_load_save_test.octbin.gz p4est_load_save_test.pvtu p4est_load_save_test.visit p4est_load_save_test_*.vtu \
                p4est_load_save_test.p4est p4est_load_save_test_2.pvtu	p4est_load_save_test_2.visit p4est_load_save_test_2*.vtu

  p4est_neighbor_iterator_test_SOURCES =  p4est_neighbor_iterator_test.cpp 
  p4est_neighbor_iterator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_neighbor_iterator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_neighbor_iterator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_neighbor_iterator_test.octbin.gz p4est_neighbor_iterator_test.pvtu p4est_neighbor_iterator_test.visit p4est_neighbor_iterator_test_*.vtu

  p4est_estimator_test_1_SOURCES =  p4est_estimator_test_1.cpp 
  p4est_estimator_test_1_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_estimator_test_1_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_estimator_test_1_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_estimator_test_1*.octbin.gz p4est_estimator_test_1*.pvtu p4est_estimator_test_1*.visit p4est_estimator_test_1*.vtu

  p4est_up_test_SOURCES =  p4est_up_test.cpp 
  p4est_up_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_up_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_up_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
  CLEANFILES += p4est_up_test.octbin.gz p4est_up_test*.pvtu p4est_up_test*.visit p4est_up_test_*.vtu

  p4est_test_discontinuous_SOURCES =  p4est_test_discontinuous.cpp
  p4est_test_discontinuous_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_test_discontinuous_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_test_discontinuous_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_test_discontinuous*.octbin.gz p4est_test_discontinuous*.pvtu p4est_test_discontinuous*.visit p4est_test_discontinuous*.vtu

  p4est_estimator_test_1_uniform_SOURCES =  p4est_estimator_test_1_uniform.cpp 
  p4est_estimator_test_1_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_estimator_test_1_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_estimator_test_1_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_estimator_test_1*.octbin.gz p4est_estimator_test_1*.pvtu p4est_estimator_test_1*.visit p4est_estimator_test_1*.vtu

  p4est_metrics_test_1_SOURCES =  p4est_metrics_test_1.cpp 
  p4est_metrics_test_1_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_metrics_test_1_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_metrics_test_1_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_metrics_test_1*.octbin.gz p4est_metrics_test_1*.pvtu p4est_metrics_test_1*.visit p4est_metrics_test_1*.vtu

  p4est_dr_test_1_marker_SOURCES =  p4est_dr_test_1_marker.cpp 
  p4est_dr_test_1_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_1_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_1_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_1_marker*.octbin.gz p4est_dr_test_1_marker*.pvtu p4est_dr_test_1_marker*.visit p4est_dr_test_1_marker*.vtu

  p4est_dr_test_1_metrics_SOURCES =  p4est_dr_test_1_metrics.cpp 
  p4est_dr_test_1_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_1_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_1_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_1_metrics*.octbin.gz p4est_dr_test_1_metrics*.pvtu p4est_dr_test_1_metrics*.visit p4est_dr_test_1_metrics*.vtu

  p4est_dr_test_2_marker_SOURCES =  p4est_dr_test_2_marker.cpp 
  p4est_dr_test_2_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_2_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_2_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_2_marker*.octbin.gz p4est_dr_test_2_marker*.pvtu p4est_dr_test_2_marker*.visit p4est_dr_test_2_marker*.vtu

  p4est_dr_test_2_metrics_SOURCES =  p4est_dr_test_2_metrics.cpp 
  p4est_dr_test_2_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_2_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_2_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_2_metrics*.octbin.gz p4est_dr_test_2_metrics*.pvtu p4est_dr_test_2_metrics*.visit p4est_dr_test_2_metrics*.vtu

  p4est_dr_test_3_marker_SOURCES =  p4est_dr_test_3_marker.cpp 
  p4est_dr_test_3_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_3_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_3_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_3_marker*.octbin.gz p4est_dr_test_3_marker*.pvtu p4est_dr_test_3_marker*.visit p4est_dr_test_3_marker*.vtu

  p4est_dr_test_3_metrics_SOURCES =  p4est_dr_test_3_metrics.cpp 
  p4est_dr_test_3_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_3_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_3_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_3_metrics*.octbin.gz p4est_dr_test_3_metrics*.pvtu p4est_dr_test_3_metrics*.visit p4est_dr_test_3_metrics*.vtu

  p4est_dr_test_4_marker_SOURCES =  p4est_dr_test_4_marker.cpp 
  p4est_dr_test_4_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_4_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_4_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_4_marker*.octbin.gz p4est_dr_test_4_marker*.pvtu p4est_dr_test_4_marker*.visit p4est_dr_test_4_marker*.vtu

  p4est_dr_test_4_metrics_SOURCES =  p4est_dr_test_4_metrics.cpp 
  p4est_dr_test_4_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_dr_test_4_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_dr_test_4_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_dr_test_4_metrics*.octbin.gz p4est_dr_test_4_metrics*.pvtu p4est_dr_test_4_metrics*.visit p4est_dr_test_4_metrics*.vtu

  p4est_adr_test_1_marker_SOURCES =  p4est_adr_test_1_marker.cpp
  p4est_adr_test_1_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_adr_test_1_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_adr_test_1_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_adr_test_1_marker*.octbin.gz p4est_adr_test_1_marker*.pvtu p4est_adr_test_1_marker*.visit p4est_adr_test_1_marker*.vtu

  p4est_adr_test_1_metrics_SOURCES =  p4est_adr_test_1_metrics.cpp
  p4est_adr_test_1_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_adr_test_1_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_adr_test_1_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_adr_test_1_metrics*.octbin.gz p4est_adr_test_1_metrics*.pvtu p4est_adr_test_1_metrics*.visit p4est_adr_test_1_metrics*.vtu

  p4est_adr_test_2_marker_SOURCES =  p4est_adr_test_2_marker.cpp 
  p4est_adr_test_2_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_adr_test_2_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_adr_test_2_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_adr_test_2_marker*.octbin.gz p4est_adr_test_2_marker*.pvtu p4est_adr_test_2_marker*.visit p4est_adr_test_2_marker*.vtu

  p4est_adr_test_2_metrics_SOURCES =  p4est_adr_test_2_metrics.cpp 
  p4est_adr_test_2_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
  p4est_adr_test_2_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_adr_test_2_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
  CLEANFILES += p4est_adr_test_2_metrics*.octbin.gz p4est_adr_test_2_metrics*.pvtu p4est_adr_test_2_metrics*.visit p4est_adr_test_2_metrics*.vtu
endif

TESTS =

if HAVE_P4EST

TESTS += p4est_quad_index_test p4est_conn_test_hanging p4est_hanging_test_2 p4est_conn_test \
         p4est_refine_segment_test p4est_hanging_test p4est_coarsen_test \
         p4est_iterator_test p4est_ref_test p4est_load_save_test p4est_edge_id_test \
         p4est_operator_test_2 p4est_operator_test_3 p4est_operator_segment_test \
         p4est_neighbor_iterator_test p4est_estimator_test_1 p4est_up_test \
         p4est_ref_lnodes_test p4est_test_discontinuous p4est_estimator_test_1_uniform \
         p4est_metrics_test_1

endif

TESTS += io_test bgs_test sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 \
	 bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc \
         bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 \
         bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc \
         bim3_timedependent_test3 bim3_timedependent_test3_osc 

if HAVE_LIS
  TESTS += linear_solver_test linear_solver_test_repeated \
           lis_stationary_test lis_timedependent_test nonlinear_solver_eq_test\
           nonlinear_solver_test nonlinear_solver_test2 nonlinear_solver_test3 nonlinear_solver_test4
endif

