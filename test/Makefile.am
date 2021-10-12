BIMMUMPSLIB = ../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../src/libbimp4est.la @P4EST_LIB@
BIMLINALGLIB = ../src/libbimlinalg.la 
BIMIOLIB    = @OCTAVE_FILE_IO_LIB@ @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../src/libbim.la @BLAS_LIB@
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@ 
AUTOMAKE_OPTIONS = nostdinc

CLEANFILES = 


check_PROGRAMS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 bim3_stationary_test linear_solver_test_simple \
                 bim3_stationary_test2 bim3_stationary_test2_osc bim3_stationary_test3 \
                 bim3_stationary_test3_osc bim3_stationary_test4 bim3_timedependent_test \
                 bim3_timedependent_test2 bim3_timedependent_test2_osc bim3_timedependent_test3 \
                 bim3_timedependent_test3_osc 

if HAVE_LIS
  check_PROGRAMS += lis_stationary_test lis_timedependent_test linear_solver_test \
                    linear_solver_test_repeated nonlinear_solver_test nonlinear_solver_test2 \
                    nonlinear_solver_test3  nonlinear_solver_eq_test #nonlinear_solver_test4
endif

if HAVE_P4EST

  check_PROGRAMS += p4est_quad_index_test p4est_conn_test_hanging p4est_hanging_test_2 \
                    p4est_conn_test p4est_hanging_test p4est_coarsen_test \
                    p4est_iterator_test  p4est_load_save_test p4est_edge_id_test \
                    p4est_neighbor_iterator_test p4est_up_test \
                    tmesh_3d_test p4est_load_save_test_3d p4est_up_test_3d \
                    p4est_refine_segment_test_3d p4est_refine_segment_test_3d_mpi.mpi p4est_ref_test_3d \
                    p4est_coarsen_test_3d p4est_hanging_test_3d p4est_operator_test_3d

endif

# schur_test_SOURCES = schur_test.cpp
# CLEANFILES += 
# schur_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB)

sparse_test_SOURCES = sparse_test.cpp
CLEANFILES += first_solve.txt second_solve.txt
sparse_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += SG.m mesh_out.m

bim3_test2_SOURCES = bim3_test2.cpp
bim3_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += SG2.m mesh_out.m

bim3_test3_SOURCES = bim3_test3.cpp
bim3_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += SG3.m mesh_out.m

bim3_test4_SOURCES = bim3_test4.cpp
bim3_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += UP.m mesh_out.m

bim3_stationary_test_SOURCES = bim3_stationary_test.cpp
bim3_stationary_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution.txt

bim3_stationary_test2_SOURCES = bim3_stationary_test2.cpp
bim3_stationary_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution2.txt mesh_out_cube.m

bim3_stationary_test2_osc_SOURCES = bim3_stationary_test2_osc.cpp
bim3_stationary_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution2_OSC.txt

bim3_stationary_test3_SOURCES = bim3_stationary_test3.cpp
bim3_stationary_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution3.txt

bim3_stationary_test3_osc_SOURCES = bim3_stationary_test3_osc.cpp
bim3_stationary_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution3_OSC.txt

bim3_stationary_test4_SOURCES = bim3_stationary_test4.cpp
bim3_stationary_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += solution4.txt

bim3_timedependent_test_SOURCES = bim3_timedependent_test.cpp
bim3_timedependent_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += Solution_TimeTest1.txt

bim3_timedependent_test2_SOURCES = bim3_timedependent_test2.cpp
bim3_timedependent_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += Solution_TimeTest2.txt

bim3_timedependent_test2_osc_SOURCES = bim3_timedependent_test2_osc.cpp
bim3_timedependent_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += Solution_TimeTest2_osc.txt

bim3_timedependent_test3_SOURCES = bim3_timedependent_test3.cpp
bim3_timedependent_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += Solution_TimeTest3.txt

bim3_timedependent_test3_osc_SOURCES = bim3_timedependent_test3_osc.cpp
bim3_timedependent_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB)
CLEANFILES += Solution_TimeTest3_osc.txt

if HAVE_LIS

  lis_stationary_test_SOURCES = lis_stationary_test.cpp
  lis_stationary_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += LIS_solution.txt

  lis_timedependent_test_SOURCES = lis_timedependent_test.cpp
  lis_timedependent_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += Lis_Solution_TimeTest.txt

  linear_solver_test_SOURCES = linear_solver_test.cpp
  linear_solver_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMLINALGLIB) $(BIMLISLIB) $(BIMMUMPSLIB) 
  CLEANFILES += MUMPS_solution.txt LIS_solution.txt

  linear_solver_test_simple_SOURCES = linear_solver_test_simple.cpp
  linear_solver_test_simple_LDADD =  $(BIMMUMPSLIB) @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
  CLEANFILES += MUMPS_solution_simple.txt LIS_solution_simple.txt

  linear_solver_test_repeated_SOURCES = linear_solver_test_repeated.cpp
  linear_solver_test_repeated_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)

  nonlinear_solver_test_SOURCES =  nonlinear_solver_test.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test_LDADD =  $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += NonLinear_Plaplacian_Test1.txt

  nonlinear_solver_test2_SOURCES =  nonlinear_solver_test2.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += NonLinear_Plaplacian_Test2.txt

  nonlinear_solver_test3_SOURCES =  nonlinear_solver_test3.cpp plaplacian_class.h plaplacian_class.cpp
  nonlinear_solver_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += NonLinear_Plaplacian_Test3.txt

  # nonlinear_solver_test4_SOURCES =  nonlinear_solver_test4.cpp plaplacian_class.h plaplacian_class.cpp
  # nonlinear_solver_test4_LDADD = $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLIB)
  # CLEANFILES += NonLinear_Plaplacian_Test4.txt

  nonlinear_solver_eq_test_SOURCES =  nonlinear_solver_eq_test.cpp equation_class.h equation_class.cpp
  nonlinear_solver_eq_test_LDADD =  $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) $(BIMLINALGLIB) $(BIMIOLIB)
  CLEANFILES += NonLinear_Equation_Test.txt

endif

if HAVE_P4EST

  p4est_quad_index_test_SOURCES =  p4est_quad_index_test.cpp 
  p4est_quad_index_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_quad_index_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_quad_index_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_quad_index_test*.octbin.gz p4est_quad_index_test.pvtu p4est_quad_index_test.visit p4est_quad_index_test_*.vtu

  p4est_conn_test_hanging_SOURCES =  p4est_conn_test_hanging.cpp 
  p4est_conn_test_hanging_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_conn_test_hanging_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_conn_test_hanging_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_conn_test_hanging*.octbin.gz p4est_conn_test_hanging.pvtu p4est_conn_test_hanging.visit p4est_conn_test_hanging_*.vtu

  p4est_hanging_test_2_SOURCES =  p4est_hanging_test_2.cpp
  p4est_hanging_test_2_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_hanging_test_2_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_hanging_test_2_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_hanging_test_2.octbin.gz p4est_hanging_test_2.pvtu p4est_hanging_test_2.visit p4est_hanging_test_2_*.vtu

  p4est_conn_test_SOURCES =  p4est_conn_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_conn_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_conn_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_conn_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_conn_test*.octbin.gz p4est_conn_test.pvtu p4est_conn_test.visit p4est_conn_test_*.vtu

  p4est_hanging_test_SOURCES =  p4est_hanging_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_hanging_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_hanging_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_hanging_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_hanging_test.octbin.gz p4est_hanging_test.pvtu p4est_hanging_test.visit p4est_hanging_test_*.vtu

  p4est_coarsen_test_SOURCES =  p4est_coarsen_test.cpp  mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_coarsen_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_coarsen_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_coarsen_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_coarsen_test.octbin.gz p4est_coarsen_test.pvtu p4est_coarsen_test.visit p4est_coarsen_test_*.vtu

  p4est_iterator_test_SOURCES =  p4est_iterator_test.cpp 
  p4est_iterator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_iterator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_iterator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_iterator_test.octbin.gz p4est_iterator_test.pvtu p4est_iterator_test.visit p4est_iterator_test_*.vtu

  p4est_edge_id_test_SOURCES =  p4est_edge_id_test.cpp 
  p4est_edge_id_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_edge_id_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_edge_id_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_edge_id_test.octbin.gz p4est_edge_id_test.pvtu p4est_edge_id_test.visit p4est_edge_id_test_*.vtu

  p4est_load_save_test_SOURCES =  p4est_load_save_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
  p4est_load_save_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_load_save_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_load_save_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_load_save_test.octbin.gz p4est_load_save_test.pvtu p4est_load_save_test.visit p4est_load_save_test_*.vtu \
                p4est_load_save_test.p4est p4est_load_save_test_2.pvtu	p4est_load_save_test_2.visit p4est_load_save_test_2*.vtu

  p4est_neighbor_iterator_test_SOURCES =  p4est_neighbor_iterator_test.cpp 
  p4est_neighbor_iterator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_neighbor_iterator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_neighbor_iterator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_neighbor_iterator_test.octbin.gz p4est_neighbor_iterator_test.pvtu p4est_neighbor_iterator_test.visit p4est_neighbor_iterator_test_*.vtu

  p4est_up_test_SOURCES =  p4est_up_test.cpp 
  p4est_up_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_up_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_up_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_up_test.octbin.gz p4est_up_test*.pvtu p4est_up_test*.visit p4est_up_test_*.vtu

  tmesh_3d_test_SOURCES =  tmesh_3d_test.cpp 
  tmesh_3d_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  tmesh_3d_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  tmesh_3d_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += tmesh_3d_out.txt

  p4est_load_save_test_3d_SOURCES =  p4est_load_save_test_3d.cpp mosfet_connectivity_3d.cpp
  p4est_load_save_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_load_save_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_load_save_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_load_save_test_3d.octbin.gz p4est_load_save_test_3d.p8est \
                p4est_load_save_test_3d*.octbin.gz \
                p4est_load_save_test_3d.pvtu p4est_load_save_test_3d.visit \
                p4est_load_save_test_3d_*.vtu p4est_load_save_test_3d_2.pvtu \
                p4est_load_save_test_3d_2.visit p4est_load_save_test_3d_2*.vtu

  p4est_up_test_3d_SOURCES =  p4est_up_test_3d.cpp 
  p4est_up_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_up_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_up_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_up_test_3d.octbin.gz p4est_up_test_3d*.pvtu \
                p4est_up_test_3d*.visit p4est_up_test_3d_*.vtu

  p4est_refine_segment_test_3d_SOURCES =  p4est_refine_segment_test_3d.cpp
  p4est_refine_segment_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_refine_segment_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_refine_segment_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_unitcube.octbin.gz p4est_refine_segment_test_3d*.pvtu \
                p4est_refine_segment_test_3d*.visit p4est_refine_segment_test_3d_*.vtu

  p4est_refine_segment_test_3d_mpi_mpi_SOURCES =  p4est_refine_segment_test_3d.cpp
  p4est_refine_segment_test_3d_mpi_mpi_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_refine_segment_test_3d_mpi_mpi_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_refine_segment_test_3d_mpi_mpi_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_unitcube.octbin.gz p4est_refine_segment_test_3d*.pvtu \
                p4est_refine_segment_test_3d*.visit p4est_refine_segment_test_3d_*.vtu

  p4est_ref_test_3d_SOURCES =  p4est_ref_test_3d.cpp mosfet_doping_3d.cpp mosfet_connectivity_3d.cpp
  p4est_ref_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_ref_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_ref_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_ref_test_3d.octbin*.gz p4est_ref_test_3d*.pvtu \
                p4est_ref_test_3d*.visit p4est_ref_test_3d_*.vtu


  p4est_coarsen_test_3d_SOURCES =  p4est_coarsen_test_3d.cpp  \
                                   mosfet_doping_3d.cpp \
                                   mosfet_connectivity_3d.cpp
  p4est_coarsen_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_coarsen_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_coarsen_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_coarsen_test_3d.octbin.gz p4est_coarsen_test_3d.pvtu \
                p4est_coarsen_test_3d.visit p4est_coarsen_test_3d_*.vtu 

  p4est_hanging_test_3d_SOURCES =  p4est_hanging_test_3d.cpp \
                                   mosfet_doping_3d.cpp \
                                   mosfet_connectivity_3d.cpp
  p4est_hanging_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_hanging_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_hanging_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_hanging_test_3d.octbin.gz p4est_hanging_test_3d.pvtu \
                p4est_hanging_test_3d.visit p4est_hanging_test_3d_*.vtu

  p4est_operator_test_3d_SOURCES =  p4est_operator_test_3d.cpp
  p4est_operator_test_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
  p4est_operator_test_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
  p4est_operator_test_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
  CLEANFILES += p4est_operator_test_3d*.octbin.gz p4est_operator_test_3d.pvtu \
                p4est_operator_test_3d.visit p4est_operator_test_3d_*.vtu

endif

TEST_EXTENSIONS = .m .sh .mpi
M_LOG_COMPILER = octave
SH_LOG_COMPILER = bash
MPIEXEC ?= mpiexec
MPI_LOG_COMPILER = $(MPIEXEC)
AM_M_LOG_FLAGS = -q --no-gui --path $(top_srcdir)/script/m:$(top_srcdir)/data
AM_MPI_LOG_FLAGS = -np 2 

TESTS =

if HAVE_P4EST

TESTS += p4est_quad_index_test p4est_conn_test_hanging p4est_hanging_test_2 p4est_conn_test \
         p4est_hanging_test p4est_coarsen_test \
         p4est_iterator_test p4est_load_save_test p4est_edge_id_test \
         p4est_neighbor_iterator_test p4est_up_test \
         tmesh_3d_test p4est_load_save_test_3d p4est_up_test_3d \
         p4est_refine_segment_test_3d p4est_refine_segment_test_3d_mpi.mpi p4est_ref_test_3d \
         p4est_coarsen_test_3d p4est_hanging_test_3d p4est_operator_test_3d linear_solver_test_simple

endif

TESTS += sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 \
	 bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc \
         bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 \
         bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc \
         bim3_timedependent_test3 bim3_timedependent_test3_osc 

if HAVE_LIS
TESTS += linear_solver_test linear_solver_test_repeated \
          lis_stationary_test lis_timedependent_test nonlinear_solver_eq_test\
          nonlinear_solver_test nonlinear_solver_test2 nonlinear_solver_test3 #nonlinear_solver_test4
endif

