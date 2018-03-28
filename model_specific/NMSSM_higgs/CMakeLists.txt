add_library(model_specific_NMSSM_higgs
  nmssm_twoloophiggs.cpp
  nmssm2loop.f)

target_include_directories(model_specific_NMSSM_higgs
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(model_specific_NMSSM_higgs
  PUBLIC FlexibleSUSY::flexisusy FlexibleSUSY::model_specific::MSSM_higgs)

add_library(FlexibleSUSY::model_specific::NMSSM_higgs
  ALIAS model_specific_NMSSM_higgs)
