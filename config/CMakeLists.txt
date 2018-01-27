string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%SZ UTC" UTC)

configure_file(config_cmake.h.in config.h)

add_library(config INTERFACE)
target_include_directories(config INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(FlexibleSUSY::config ALIAS config)
