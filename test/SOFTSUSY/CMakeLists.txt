set(WITH_SOFTSUSY 1 CACHE INTERNAL WITH_SOFTSUSY)

set(SOFTSUSY_SRC
  conversion.cpp
  def.cpp
  diagonalization.cpp
  flavoursoft.cpp
  linalg.cpp
  lowe_legacy.cpp
  mssmUtils.cpp
  nmssmUtils.cpp
  nmssmsoftpars.cpp
  nmssmsoftsusy.cpp
  nmssmsusy.cpp
  nmssm1loop.f
  numerics_legacy.cpp
  physpars.cpp
  rge.cpp
  rk_legacy.cpp
  susy.cpp
  tensor.cpp
  utils.cpp
)

if(ENABLE_TWO_SCALE_SOLVER)
  set(SOFTSUSY_SRC
    ${SOFTSUSY_SRC}
    SoftsusyMSSM_two_scale.cpp
    SoftsusyMSSM_two_scale_convergence_tester.cpp
    SoftsusyMSSM_two_scale_initial_guesser.cpp
    SoftsusyMSSM_two_scale_susy_scale_constraint.cpp
    SoftsusyMSSM_two_scale_low_scale_constraint.cpp
    SoftsusyMSSM_two_scale_sugra_constraint.cpp
    SoftsusyNMSSM_two_scale.cpp
    SoftsusyNMSSM_two_scale_convergence_tester.cpp
    SoftsusyNMSSM_two_scale_initial_guesser.cpp
    SoftsusyNMSSM_two_scale_low_scale_constraint.cpp
    SoftsusyNMSSM_two_scale_sugra_constraint.cpp
    SoftsusyNMSSM_two_scale_susy_scale_constraint.cpp
  )
endif()

add_library(SOFTSUSY ${SOFTSUSY_SRC})

target_include_directories(SOFTSUSY
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SOFTSUSY
  PUBLIC
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::MSSM_higgs
  FlexibleSUSY::model_specific::NMSSM_higgs
)

add_library(FlexibleSUSY::SOFTSUSY ALIAS SOFTSUSY)

add_executable(run_softpoint.x run_softpoint.cpp)
target_link_libraries(run_softpoint.x FlexibleSUSY::SOFTSUSY)

add_executable(run_softsusy.x run_softsusy.cpp)
target_link_libraries(run_softsusy.x FlexibleSUSY::SOFTSUSY)
