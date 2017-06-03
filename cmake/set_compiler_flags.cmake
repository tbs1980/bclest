include(CheckCXXCompilerFlag)

macro(add_cxx_compiler_flag FLAG ARG_BUILD_TYPE)
  string(REGEX REPLACE "-" "" SFLAG1 ${FLAG})
  string(REGEX REPLACE "\\+" "p" SFLAG ${SFLAG1})
  check_cxx_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})

  if(COMPILER_SUPPORT_${SFLAG})
    if(${ARG_BUILD_TYPE} MATCHES "DEBUG")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${FLAG}")
    elseif(${ARG_BUILD_TYPE} MATCHES "RELEASE")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${FLAG}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
  endif()
endmacro()

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_cxx_compiler_flag("-g3" "DEBUG")
    add_cxx_compiler_flag("-g0" "RELEASE")
    add_cxx_compiler_flag("-O3" "RELEASE")

    add_cxx_compiler_flag("-std=c++11" "")
    add_cxx_compiler_flag("-pedantic" "")
    add_cxx_compiler_flag("-Wall" "")
    add_cxx_compiler_flag("-Wextra" "")
  endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_cxx_compiler_flag("-g3" "DEBUG")
    
    # for code coverage
    add_cxx_compiler_flag("-O0" "DEBUG")
    add_cxx_compiler_flag("-fprofile-arcs" "DEBUG")
    add_cxx_compiler_flag("-ftest-coverage" "DEBUG")

    add_cxx_compiler_flag("-g0" "RELEASE")
    add_cxx_compiler_flag("-O3" "RELEASE")

    add_cxx_compiler_flag("-std=c++11" "")
    add_cxx_compiler_flag("-pedantic" "")
    add_cxx_compiler_flag("-Wall" "")
    add_cxx_compiler_flag("-Wextra" "")
  endif()
endif()