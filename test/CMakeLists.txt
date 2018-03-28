add_library(flexiblesusy_test
  error_count.cpp
  run_cmd.cpp
  stopwatch.cpp)

target_include_directories(flexiblesusy_test
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(flexiblesusy_test
  PUBLIC FlexibleSUSY::flexisusy ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_library(FlexibleSUSY::flexiblesusy_test ALIAS flexiblesusy_test)

function(add_flexiblesusy_bin_test name)
  if(ENABLE_COMPILE)
    add_executable(
      ${name}
      ${name}.cpp
      )
    target_link_libraries(
      ${name}
      FlexibleSUSY::flexiblesusy_test
      FlexibleSUSY::model_specific::MSSM_higgs
      FlexibleSUSY::model_specific::MSSM_thresholds
      FlexibleSUSY::model_specific::NMSSM_higgs
      FlexibleSUSY::model_specific::SM
      FlexibleSUSY::model_specific::SplitMSSM
      )
    add_test(NAME ${name} COMMAND ${name})
  endif()
endfunction()

add_flexiblesusy_bin_test(test_array_view)
add_flexiblesusy_bin_test(test_cast_model)
add_flexiblesusy_bin_test(test_ckm)
add_flexiblesusy_bin_test(test_logger)
add_flexiblesusy_bin_test(test_derivative)
add_flexiblesusy_bin_test(test_effective_couplings)
add_flexiblesusy_bin_test(test_eigen_utils)
add_flexiblesusy_bin_test(test_ewsb_solver)
add_flexiblesusy_bin_test(test_fixed_point_iterator)
add_flexiblesusy_bin_test(test_goldstones)
add_flexiblesusy_bin_test(test_gsl_vector)
add_flexiblesusy_bin_test(test_linalg2)
add_flexiblesusy_bin_test(test_minimizer)
add_flexiblesusy_bin_test(test_namespace_collisions)
add_flexiblesusy_bin_test(test_mssm_twoloop_mb)
add_flexiblesusy_bin_test(test_mssm_twoloop_mt)
add_flexiblesusy_bin_test(test_MSSM_2L_limits)
add_flexiblesusy_bin_test(test_numerics)
add_flexiblesusy_bin_test(test_problems)
add_flexiblesusy_bin_test(test_pv)
add_flexiblesusy_bin_test(test_raii)
add_flexiblesusy_bin_test(test_root_finder)
add_flexiblesusy_bin_test(test_scan)
add_flexiblesusy_bin_test(test_sminput)
add_flexiblesusy_bin_test(test_slha_io)
add_flexiblesusy_bin_test(test_sum)
add_flexiblesusy_bin_test(test_threshold_corrections)
add_flexiblesusy_bin_test(test_threshold_loop_functions)
add_flexiblesusy_bin_test(test_spectrum_generator_settings)
add_flexiblesusy_bin_test(test_which)
add_flexiblesusy_bin_test(test_wrappers)

if(WITH_GM2Calc)
  add_flexiblesusy_bin_test(test_gm2calc)
  add_flexiblesusy_bin_test(test_MSSMNoFV_onshell)
endif()
