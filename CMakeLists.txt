project(Eigenvalues)
cmake_minimum_required(VERSION 3.0)
if(ENABLE_COVERAGE)
    include(Coverage)
endif()

# boost
message(STATUS "Checking for Boost...")
set(Boost_ADDITIONAL_VERSIONS "1.62;1.61;1.60;1.59;1.57;1.56;1.55;1.41.0;1.41;1.42;1.43;1.44;1.45;1.46;1.46.1;1.47;1.50;1.51;1.54")
if(ENABLE_STATIC_EXTERNALS)
    set(Boost_USE_STATIC_LIBS TRUE)
endif()
set(Boost_DETAILED_FAILURE_MSG TRUE)
if(NOT BOOST_ROOT AND NOT DEFINED $ENV{BOOST_ROOT})
    if(NUTO_EXTERNAL_LIBRARIES_ROOT)
        set(BOOST_ROOT ${NUTO_EXTERNAL_LIBRARIES_ROOT})
    endif()
endif()

find_package(Boost 1.5.8.9 COMPONENTS unit_test_framework REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost version = ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "BOOST_INCLUDE_DIR = ${BOOST_INCLUDE_DIR}")
message(STATUS "BOOST_INCLUDE_DIRS = ${BOOST_INCLUDE_DIRS}")

#test Environment
include(CTest)

add_executable(TestEigenvalue TestEigenvalue.cpp)
target_link_libraries(TestEigenvalue ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(Test::Eigenvalue TestEigenvalue --log_level=message)
