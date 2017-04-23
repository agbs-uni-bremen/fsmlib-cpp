#set the minimum cmake version required
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.11)

#set the name of the project
PROJECT (fsm)

#####################################################################
#enable warning for MSVC and gcc

if (MSVC)
	# Force to always compile with W4
	if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string (REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
	endif ()
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
add_subdirectory (example)
add_subdirectory (generator)

if(gui)
	add_subdirectory (window)

#copy every fsm file to the build directory
configure_file (../resources/fsm.fsm									window/fsm.fsm									COPYONLY)
configure_file (../resources/fsm.in										window/fsm.in										COPYONLY)
configure_file (../resources/fsm.out									window/fsm.out									COPYONLY)
configure_file (../resources/fsm.state									window/fsm.state									COPYONLY)
configure_file (../resources/fsma.fsm									window/fsma.fsm									COPYONLY)
configure_file (../resources/fsma.in									window/fsma.in									COPYONLY)
configure_file (../resources/fsma.out									window/fsma.out									COPYONLY)
configure_file (../resources/fsma.state									window/fsma.state									COPYONLY)
configure_file (../resources/fsmb.fsm									window/fsmb.fsm									COPYONLY)
configure_file (../resources/fsmb.in									window/fsmb.in									COPYONLY)
configure_file (../resources/fsmb.out									window/fsmb.out									COPYONLY)
configure_file (../resources/fsmb.state									window/fsmb.state									COPYONLY)
configure_file (../resources/fsmGillA7.fsm								window/fsmGillA7.fsm								COPYONLY)
configure_file (../resources/fsmGillA7.in								window/fsmGillA7.in								COPYONLY)
configure_file (../resources/fsmGillA7.out								window/fsmGillA7.out								COPYONLY)
configure_file (../resources/fsmGillA7.state							window/fsmGillA7.state							COPYONLY)
endif(gui)
