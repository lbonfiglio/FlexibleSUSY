string(REPLACE "," ";" MODELS "${WITH_MODELS}")

set(CMakeLists_template
  "\
# Generated CMakeLists.txt file for the model @MODEL@.
# Generated at @DATE_AND_TIME@.

if(ENABLE_META)
  Mathematica_ADD_CUSTOM_TARGET(
    # OUTPUT 00_DELETE_ME_TO_RERUN_METACODE
    generate-@MODEL@-src
    SCRIPT start.m
    DEPENDS
      \"start.m\"
      \"FlexibleSUSY.m\"
      \\\${META_SRC}
      \\\${TEMPLATE_SRC}
      @SARAH_MODEL_FILES@
    KERNEL_OPTIONS -noprompt
    COMMENT \"Generating C++ source code for @MODEL@\"
  )
endif()

if(ENABLE_COMPILE)
  file(GLOB @MODEL@_SRC *.cpp)

  add_library(@MODEL@ \\\${@MODEL@_SRC})
  target_include_directories(@MODEL@ PUBLIC @DIR@)
  target_link_libraries(@MODEL@ PUBLIC FlexibleSUSY::flexisusy)
  add_library(FlexibleSUSY::@MODEL@ ALIAS @MODEL@)

  add_executable(run_@MODEL@ run_@MODEL@.cpp)
  target_include_directories(run_@MODEL@ PRIVATE @DIR@)
  target_link_libraries(run_@MODEL@ PRIVATE FlexibleSUSY::@MODEL@)

  if(ENABLE_LIBRARYLINK AND Mathematica_FOUND)
    Mathematica_ADD_LIBRARY(@MODEL@_libraylink @MODEL@_libraylink.cpp)
  endif()
endif()
"
)

macro(fill_template _content _template _model)
  string(TIMESTAMP DATE "%a %d %b %Y %H:%M:%SZ UTC" UTC)
  set(${_content} "${_template}")
  string(REPLACE "@MODEL@" "${_model}" ${_content} "${${_content}}")
  string(REPLACE "@DATE_AND_TIME@" "${DATE}" ${_content} "${${_content}}")
  string(REPLACE "@DIR@" "\${CMAKE_CURRENT_SOURCE_DIR}" ${_content} "${${_content}}")
  # TODO: fill SARAH_MODEL_FILES
  string(REPLACE "@SARAH_MODEL_FILES@" "${SARAH_MODEL_FILES}" ${_content} "${${_content}}")
endmacro()

function(update_cmake_file model)
  # check whether the CMakeLists.txt file should be updated
  # TODO: set MUST_UPDATE depending on meta code / templates have changed
  set(MUST_UPDATE ON)

  if(MUST_UPDATE OR FORCE_RECREATE)
    # create CMakeLists.txt
    fill_template(content "${CMakeLists_template}" "${model}")
    file(WRITE "${model}/CMakeLists.txt" "${content}")
    # TODO: create start.m
    # TODO: copy FlexibleSUSY.m from model_files/
    # TODO: copy LesHouches.* from model_files/
  endif()
endfunction()

function(process_model model)
  message(STATUS "Processing model: ${model}")
  update_cmake_file(${model})
  add_subdirectory(${model})
endfunction()

foreach(model ${MODELS})
  process_model(${model})
endforeach()
