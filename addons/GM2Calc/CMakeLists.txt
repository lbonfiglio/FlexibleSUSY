if(ENABLE_GM2Calc)
  set(WITH_GM2Calc 1 CACHE INTERNAL WITH_GM2Calc)

  add_library(GM2Calc
    ffunctions.cpp
    gm2_1loop.cpp
    gm2_2loop.cpp
    gm2_mb.cpp
    gm2_slha_io.cpp
    gm2_uncertainty.cpp
    gm2calc_interface.cpp
    MSSMNoFV_onshell.cpp
    MSSMNoFV_onshell_mass_eigenstates.cpp
    MSSMNoFV_onshell_physical.cpp
    MSSMNoFV_onshell_problems.cpp
    MSSMNoFV_onshell_soft_parameters.cpp
    MSSMNoFV_onshell_susy_parameters.cpp)

  target_link_libraries(GM2Calc PUBLIC FlexibleSUSY::flexisusy)
  target_include_directories(GM2Calc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  add_library(FlexibleSUSY::GM2Calc ALIAS GM2Calc)

  add_executable(gm2calc gm2calc.cpp)
  target_link_libraries(gm2calc GM2Calc)
  target_include_directories(gm2calc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
else()
  add_library(GM2Calc INTERFACE)
  add_library(FlexibleSUSY::GM2Calc ALIAS GM2Calc) 
  target_include_directories(GM2Calc INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
