add_library(model_specific_MSSM_higgs
  mssm_twoloophiggs.cpp
  mssm_twoloophiggs_impl.f)

target_include_directories(model_specific_MSSM_higgs
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(model_specific_MSSM_higgs
  PUBLIC FlexibleSUSY::flexisusy)

add_library(FlexibleSUSY::model_specific::MSSM_higgs
  ALIAS model_specific_MSSM_higgs)
