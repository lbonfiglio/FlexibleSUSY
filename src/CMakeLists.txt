set(FLEXISUSY_SRC
  betafunction.cpp
  build_info.cpp
  bvp_solver_problems.cpp
  ckm.cpp
  command_line_options.cpp
  composite_convergence_tester.cpp
  database.cpp
  dilog.cpp
  dilogc.f
  effective_couplings.cpp
  global_thread_pool.cpp
  gsl_utils.cpp
  gsl_vector.cpp
  lowe.cpp
  sfermions.cpp
  numerics.cpp
  numerics2.cpp
  physical_input.cpp
  pmns.cpp
  problems.cpp
  pv.cpp
  rkf_integrator.cpp
  scan.cpp
  slha_io.cpp
  spectrum_generator_problems.cpp
  spectrum_generator_settings.cpp
  string_utils.cpp
  threshold_corrections.cpp
  threshold_loop_functions.cpp
  weinberg_angle.cpp
  wrappers.cpp)

if(ENABLE_TWO_SCALE_SOLVER)
  set(FLEXISUSY_SRC
    ${FLEXISUSY_SRC}
    two_scale_running_precision.cpp
    two_scale_solver.cpp)
endif()

if(ENABLE_SEMI_ANALYTIC_SOLVER)
  set(FLEXISUSY_SRC
    ${FLEXISUSY_SRC}
    semi_analytic_solver.cpp
    two_scale_running_precision.cpp
    two_scale_solver.cpp)
endif()

add_library(flexisusy ${FLEXISUSY_SRC})

target_include_directories(flexisusy
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(flexisusy
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  GSL::gsl
  ${LAPACK_LIBRARIES}
  Threads::Threads
)

add_library(FlexibleSUSY::flexisusy ALIAS flexisusy)
