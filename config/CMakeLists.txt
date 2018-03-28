string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%SZ UTC" UTC)

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_COMMIT)

configure_file(config_cmake.h.in config.h @ONLY)

add_library(config INTERFACE)
target_include_directories(config INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")

add_library(FlexibleSUSY::config ALIAS config)
