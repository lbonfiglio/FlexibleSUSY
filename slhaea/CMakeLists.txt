add_library(slhaea INTERFACE)
target_include_directories(slhaea INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(FlexibleSUSY::slhaea ALIAS slhaea)
