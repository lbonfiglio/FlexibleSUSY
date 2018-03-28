add_library(model_specific_MSSM_thresholds
  mssm_twoloop_as.cpp
  mssm_twoloop_mb.cpp
  mssm_twoloop_mt.cpp
  mssm_twoloop_mtau.cpp)

target_include_directories(model_specific_MSSM_thresholds
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(model_specific_MSSM_thresholds
  PUBLIC FlexibleSUSY::flexisusy)

add_library(FlexibleSUSY::model_specific::MSSM_thresholds
  ALIAS model_specific_MSSM_thresholds)
