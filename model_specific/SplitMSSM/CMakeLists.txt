add_library(model_specific_SplitMSSM
  splitmssm_threeloophiggs.cpp
  splitmssm_thresholds.cpp)

target_include_directories(model_specific_SplitMSSM
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(model_specific_SplitMSSM
  PUBLIC FlexibleSUSY::flexisusy)

add_library(FlexibleSUSY::model_specific::SplitMSSM
  ALIAS model_specific_SplitMSSM)
