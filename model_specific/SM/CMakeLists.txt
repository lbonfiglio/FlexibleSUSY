add_library(model_specific_SM
  sm_fourloophiggs.cpp
  sm_threeloop_as.cpp
  sm_threeloophiggs.cpp
  sm_twoloophiggs.cpp
  standard_model.cpp
  standard_model_effective_couplings.cpp
  standard_model_physical.cpp
  standard_model_two_scale_convergence_tester.cpp
  standard_model_two_scale_low_scale_constraint.cpp
  standard_model_two_scale_model.cpp)

target_include_directories(model_specific_SM
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(model_specific_SM
  PUBLIC FlexibleSUSY::flexisusy)

add_library(FlexibleSUSY::model_specific::SM
  ALIAS model_specific_SM)
