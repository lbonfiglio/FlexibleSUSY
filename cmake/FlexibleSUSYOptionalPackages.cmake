function(fs_package_status_message msg)
  message(STATUS "${Blue}${msg}${ColourReset}")
endfunction()

macro(find_optional_package name enable_name)
  if(${enable_name})
    find_package(${name})
    if(${name}_FOUND)
      fs_package_status_message("Enabling use of ${name}")
    else()
      unset(${enable_name})
      fs_package_status_message("Disabling use of ${name}")
    endif()
  else()
    fs_package_status_message("Disabling use of ${name}")
  endif()
endmacro()

function(find_threads)
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  set(THREADS_PREFER_PTHREAD_FLAG TRUE)
  find_optional_package(Threads ENABLE_THREADS)
  if(NOT Threads::Threads)
    add_library(fs_threads INTERFACE)
    add_library(Threads::Threads ALIAS fs_threads)
  endif()
endfunction()

function(find_lapack)
  find_optional_package(LAPACK ENABLE_LAPACK)
endfunction()

function(find_odeint)
  if(ENABLE_ODEINT)
    find_file(BOOST_FS_ODEINT_1 boost/numeric/odeint.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_2 boost/numeric/odeint/algebra/algebra_dispatcher.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_3 boost/numeric/odeint/algebra/vector_space_algebra.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_4 boost/numeric/odeint/external/eigen/eigen_algebra_dispatcher.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_5 boost/numeric/odeint/external/eigen/eigen_resize.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_6 boost/numeric/odeint/integrate/integrate_adaptive.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_7 boost/numeric/odeint/stepper/generation.hpp ${Boost_INCLUDE_DIRS})
    find_file(BOOST_FS_ODEINT_8 boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp ${Boost_INCLUDE_DIRS})
    if(BOOST_FS_ODEINT_1 AND
        BOOST_FS_ODEINT_2 AND
        BOOST_FS_ODEINT_3 AND
        BOOST_FS_ODEINT_4 AND
        BOOST_FS_ODEINT_5 AND
        BOOST_FS_ODEINT_6 AND
        BOOST_FS_ODEINT_7 AND
        BOOST_FS_ODEINT_8)
      fs_package_status_message("Enabling use of Boost's odeint")
    else()
      unset(ENABLE_ODEINT)
      fs_package_status_message("Disabling use of Boost's odeint")
    endif()
  else()
    fs_package_status_message("Disabling use of Boost's odeint")
  endif()
endfunction()
