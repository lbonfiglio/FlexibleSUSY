if(ENABLE_FFLITE)
  add_library(fflite
    BcoeffAD.F
    ffca0.F
    ffcb0.F
    ffcb1.F
    ffcb2p.F
    ffcc0.F
    ffcli2.F
    ffinit.F
    ffxa0.F
    ffxb0.F
    ffxb1.F
    ffxb2p.F
    ffxli2.F
    ini.F)
else()
  add_library(fflite INTERFACE)
endif()

target_include_directories(fflite PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(FlexibleSUSY::fflite ALIAS fflite)
