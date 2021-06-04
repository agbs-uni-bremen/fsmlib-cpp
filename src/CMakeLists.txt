#set the minimum cmake version required
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.11)

#set the name of the project
PROJECT (fsm)
#####################################################################
#Setting debug macro if this is a DEBUG build

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_DEBUG_MACRO")

#####################################################################
# Adding flag, in order to auto configure the Easylogging++ logger
add_definitions(-DAUTO_INITIALIZE_EASYLOGGINGPP)

#####################################################################
# Configure logger

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DELPP_FEATURE_PERFORMANCE_TRACKING)
    add_definitions(-DELPP_LOGGING_FLAGS_FROM_ARG)
    #add_definitions(-DELPP_PERFORMANCE_MICROSECONDS)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DELPP_FEATURE_PERFORMANCE_TRACKING)
    #add_definitions(-DELPP_PERFORMANCE_MICROSECONDS)
endif()

#####################################################################
#enable warning for MSVC and gcc

if (MSVC)
	# Force to always compile with W4
	if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string (REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
	endif ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /permissive-")
elseif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

#####################################################################
#enable c++ 11 for MSVC and gcc

include (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG ("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG ("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else ()
	message (STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif ()

#####################################################################
OPTION( gui "Build with gui support" OFF)


if(gui)
	find_package (Qt5Widgets REQUIRED)
endif(gui)

#set the root source diectory as include directory
include_directories (${CMAKE_SOURCE_DIR})
include_directories (${CMAKE_SOURCE_DIR}/externals/jsoncpp-0.10.0/include)
include_directories (${CMAKE_BINARY_DIR})

#add every subdirectory to cmake
add_subdirectory (externals)
add_subdirectory (fsm)
add_subdirectory (interface)
add_subdirectory (main)
add_subdirectory (sets)
add_subdirectory (trees)
add_subdirectory (harness)
add_subdirectory (harness_adaptive)
add_subdirectory (example)
add_subdirectory (generator)
add_subdirectory (checker)
add_subdirectory (cloneable)
add_subdirectory (utils)

include(cppcheck.cmake)
