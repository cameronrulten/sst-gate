## ==============================================================================
##
##  SST-GATE simulation software CMake build system
##  Copyright (c) 2017, Cameron Rulten, University of Minnesota (2017)
##
##  CMake file for the SST-GATE optical perfermance analysis simulation
##
##  Details: This is the top-level CMakeLists.txt file
##
## ==============================================================================


## Set project name 
project(SSTGATE CXX)

## Set minimum cmake version required for software to install properly
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

## Set the software version number
set (SSTGATE_VERSION_MAJOR 1)
set (SSTGATE_VERSION_MINOR 0)
set (SSTGATE_PATCH_VERSION 0)

## Set the C++ standard
set(CMAKE_CXX_STANDARD 11)

#if CMAKE_INSTALL_PREFIX set to default (i.e. /usr/local) then change to binary directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE)
endif()

set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${SSTGATE_VERSION}" CACHE PATH "SSTGATE installation directory." FORCE)

#define the SSTGATE install directory as CMAKE_INSTALL_PREFIX/SSTGATE_VERSION
set(SSTGATE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" CACHE PATH "SSTGATE installation directory.")

#define the SSTGATE install lib, include and bin directories
set(SSTGATE_INSTALL_LIBRARY_DIR "${SSTGATE_INSTALL_DIR}/lib" CACHE PATH "SSTGATE installation directory for libraries")
set(SSTGATE_INSTALL_INCLUDE_DIR "${SSTGATE_INSTALL_DIR}/include" CACHE PATH "SSTGATE installation directory for header files")
set(SSTGATE_INSTALL_BINARY_DIR "${SSTGATE_INSTALL_DIR}/bin" CACHE PATH "SSTGATE installation directory for executables")

#save SSTGATE version information to text file in installatation directory
#set(version_outfile "${CMAKE_INSTALL_PREFIX}/${SSTGATE_VERSION}/SSTGATE-version")
set(version_outfile "${CMAKE_BINARY_DIR}/SSTGATE-version")
file(WRITE ${version_outfile} ${SSTGATE_VERSION})

## ==============================================================================
##
##  Set compiler specific definitions
##
## ==============================================================================

message(STATUS "${CMAKE_CXX_COMPILER_ID}")

# Check if the C++ compiler is C++11 compatible
#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR "Linux")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_ID MATCHES "Clang")
    # using regular Clang or AppleClang
    message(STATUS "Found Clang compiler collection")
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" COMPILER_SUPPORTS_C99)
    if(COMPILER_SUPPORTS_C99)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl -undefined dynamic_lookup")
      set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -m64")
    endif()
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")#" -Qunused-arguments -Wall -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter")#-O2 -DNDEBUG
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl -undefined dynamic_lookup") 
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -m64")
    elseif(COMPILER_SUPPORTS_CXX11 AND COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")#" -Wall ")#-O2 -DNDEBUG ##" -I/usr/local/include"
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++11 -g -Wl -s -undefined dynamic_lookup -03")
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -m64")
    else()
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" AND "${CMAKE_C_COMPILER_ID}" MATCHES "GNU")#STREQUAL
    # using GCC and G++
    message(STATUS "Found GNU compiler collection")
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" COMPILER_SUPPORTS_C99)#-std=c99
    if(COMPILER_SUPPORTS_C99)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -W -Wall -fsigned-char -fno-common")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl -undefined dynamic_lookup")
      #set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -Wl,-s")
      set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -bind_at_load -m64")
    endif()
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")# -stdlib=libc++ -Qunused-arguments -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter -fno-common")#-O2 -DNDEBUG-stdlib=libstdc++
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-s") 
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -Wl -undefined dynamic_lookup")
      set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")
    elseif(COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")# -stdlib=libc++ -Qunused-arguments -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter -fno-common")#-O2 -DNDEBUG
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-s")
      set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -Wl -undefined dynamic_lookup")
      set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")
    else()
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  endif()
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_ID MATCHES "Clang")
    # using regular Clang or AppleClang
    message(STATUS "Found Clang compiler collection")
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" COMPILER_SUPPORTS_C99)
    if(COMPILER_SUPPORTS_C99)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    endif()
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
    #CHECK_CXX_COMPILER_FLAG("-stdlib=stdlibc++" COMPILER_SUPPORTS_STDLIBCXX)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      #if(COMPILER_SUPPORTS_LIBCXX)
      # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
      #elseif(COMPILER_SUPPORTS_STDLIBCXX)
      # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=stdlibc++")
      #endif()
      #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")#" -Qunused-arguments -Wall -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter")#-O2 -DNDEBUG
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++11 -g -Wl -s -undefined dynamic_lookup") 
      #set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -m64")
    elseif(COMPILER_SUPPORTS_CXX11 AND COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")# -stdlib=libc++")#" -Wall ")#-O2 -DNDEBUG ##" -I/usr/local/include"
      #if(COMPILER_SUPPORTS_LIBCXX)
      # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
      #endif()
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -std=c++11 -g -Wl -s -undefined dynamic_lookup -03")
      #set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -m64")
    else()
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" AND "${CMAKE_C_COMPILER_ID}" MATCHES "GNU")#STREQUAL
    # using GCC and G++
    message(STATUS "Found GNU compiler collection")
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-std=c99" COMPILER_SUPPORTS_C99)#-std=c99
    if(COMPILER_SUPPORTS_C99)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -W -Wall -fsigned-char -fno-common")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl -undefined dynamic_lookup")
      #set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -Wl,-s")
      #set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -bind_at_load -m64")
    endif()
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")# -stdlib=libc++ -Qunused-arguments -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter -fno-common")#-O2 -DNDEBUG-stdlib=libstdc++
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-s") 
      #set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -Wl -undefined dynamic_lookup")
      #set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")
    elseif(COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")# -stdlib=libc++ -Qunused-arguments -Wextra -Wno-error=unused-command-line-argument -Wno-unused-parameter -fno-common")#-O2 -DNDEBUG
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
      #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-s")
      #set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -Wl -undefined dynamic_lookup")
      #set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -bind_at_load -m64")
    else()
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  endif()
endif()




## Check if the C++ compiler is C++11 compatible
#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR "Linux")
#  include(CheckCXXCompilerFlag)
#  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#  if(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#  elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#  else()
#    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#  endif()
#endif()

## Set compiler specific warnings
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g") ##" -I/usr/local/include")
#  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib")
#  set(warnings "-g -Qunused-arguments -Wall -Wextra -Wno-error=unused-command-line-argument")
#endif()
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#  set(warnings "-g -Wall -Wextra -Werror")
#endif()
#if (NOT CONFIGURED_ONCE)
#  set(CMAKE_CXX_FLAGS "${warnings}"
#   CACHE STRING "Flags used by the compiler during all build types." FORCE)
#  set(CMAKE_C_FLAGS   "${warnings}"
#    CACHE STRING "Flags used by the compiler during all build types." FORCE)
#endif()

# ## ==============================================================================
# ##
# ##  Set options for building libraries i.e. RPATHs and Archive warnings
# ##  Libraries will be built in the ${CMAKE_INSTALL_PREFIX}/lib directory
# ##
# ## ==============================================================================

# #set (CMAKE_SHARED_LINKER_FLAGS "-flat_namespace")
SET(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")


#In many cases we want to make sure that the required libraries are always found
#independently from LD_LIBRARY_PATH and the install location. To always use RPATH
#we set the following options:

#To enable RPATH usage on OSX and osMac
set(CMAKE_MACOSX_RPATH 1)
#set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${SSTGATE_VERSION}/lib")
SET(CMAKE_INSTALL_RPATH "${SSTGATE_INSTALL_LIBRARY_DIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
#LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${SSTGATE_VERSION}/lib" isSystemDir)
#LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${SSTGATE_INSTALL_DIR}/lib" isSystemDir)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${SSTGATE_INSTALL_LIBRARY_DIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   #SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${SSTGATE_VERSION}/lib")
   #SET(CMAKE_INSTALL_RPATH "${SSTGATE_INSTALL_DIR}/lib")
   SET(CMAKE_INSTALL_RPATH "${SSTGATE_INSTALL_LIBRARY_DIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

## ==============================================================================
##
##  System inspection for dependencies
##
## ==============================================================================


## ==============================================================================
##  CERN ROOT ANALYSIS FRAMEWORK 
## ==============================================================================


### Search environment variables to see if the ROOTSYS path has been set which is required!!!
#set(root_path "$ENV{ROOTSYS}")
##message("ROOT found! ${root_path}")
## You need to tell CMake where to find the ROOT installation
##   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${root_path}/etc/cmake/")
#set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${root_path}/etc/cmake/")

#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")

# Locate the ROOT package which is required and define required components and set extra compiler flags and include paths
find_package(ROOT 6.06 REQUIRED COMPONENTS
  Eve
  EG
  TreePlayer
  Geom
  Ged
  RGL
  Gui
  Core
  RIO
  Net
  Hist
  Graf
  Graf3d
  Gpad
  Tree
  Rint
  Postscript
  Matrix
  Physics
  MathCore
  MathMore
  Thread
  MultiProc
  XMLIO
  Minuit
  TMVA
  )
if(ROOT_FOUND)
	message("ROOT ${ROOT_VERSION} found")
	include_directories(${ROOT_INCLUDE_DIRS}) #add the root include directory to the build path
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}") #add the root cflags to the compiler
  #set(ROOT_DEFINITIONS "${ROOT_DEFINITTIONS} HAVE_ROOT=1") #add precompile definitions
  include(${ROOT_USE_FILE})
  message("ROOT libs: ${ROOT_LIBRARIES}")
  #list(APPEND CMAKE_MODULE_PATH "${root_path}/cmake")
  #list(APPEND CMAKE_MODULE_PATH "${root_path}/cmake/modules")
else(NOT ROOT_FOUND)
  message("You need CERN ROOT >= v6.06 for this program to work!")
endif(ROOT_FOUND)

# ## ==============================================================================
# ##  VEGAS ANALYSIS SOFTWARE
# ## ==============================================================================

# find_package(VEGAS 2.05.06.1 EXACT CONFIG REQUIRED)
# if(VEGAS_FOUND)
#   message(STATUS "VEGAS ${VEGAS_VERSION} found.")
#   include_directories(${VEGAS_INCLUDE_DIRS})
#   include(${VEGAS_LIBRARIES}/cmake/VEGASTargets.cmake)
#   message("VEGAS header files: ${VEGAS_INCLUDE_DIRS}")
#   message("VEGAS libs: ${VEGAS_LIBRARIES}")
#   #message("VEGAS cflags: ${VEGAS_CFLAGS} ${VEGAS_CFLAGS_OTHER}")
#   #message("VEGAS ldflags: ${VEGAS_LDFLAGS} ${VEGAS_LDFLAGS_OTHER}")
# else(NOT VEGAS_FOUND)
#   message(FATAL_ERROR "VEGAS not found. You need VEGAS for this program to work!")
# endif(VEGAS_FOUND)


# ### Search environment variables to see if the VEGAS path has been set which is required!!!
# #set(vegas_path "$ENV{VEGAS}")
# ###message("VEGAS found! ${vegas_path}")
# #find_package(PkgConfig REQUIRED)
# #if(PKG_CONFIG_FOUND)
# #	message(STATUS "PkgConfig found.")
# #else(NOT PKG_CONFIG_FOUND)
# #	message(STATUS "PkgConfig not found")
# #endif(PKG_CONFIG_FOUND)

# ###find_package(VEGAS 2.5.5 REQUIRED)
# #pkg_search_module(VEGAS 2.5.5 REQUIRED vegas)
# ###find_library(VEGAS SP24sharedLite ${vegas_path}/common/lib)
# #if(VEGAS_FOUND)
# #	message(STATUS "VEGAS ${PKG_CONFIG_VERSION_STRING} found! ${vegas_path}")
# #	#include_directories(${VEGAS_INCLUDE_DIRS})
# #	add_definitions(${VEGAS_CFLAGS} ${VEGAS_CFLAGS_OTHER})
# #	#set(VEGAS_DEFINITIONS "${VEGAS_DEFINITTIONS} -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -DNOHDF5 -D_OAWG -DHAVE_ROOT_6")
# #	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VEGAS_CXX_FLAGS}") #add the vegas cflags to the compiler
# #	message("VEGAS dirs: ${VEGAS_INCLUDE_DIRS} and ${VEGAS_LIBRARY_DIRS}")
# #	message("VEGAS libs: ${VEGAS_LIBRARIES} and ${VEGAS_LIBS}")
# #	message("VEGAS cflags: ${VEGAS_CFLAGS} ${VEGAS_CFLAGS_OTHER}")
# #	message("VEGAS ldflags: ${VEGAS_LDFLAGS} ${VEGAS_LDFLAGS_OTHER}")
# #else(NOT VEGAS_FOUND)
# #	message(STATUS "VEGAS not found. You need VEGAS for this program to work!")
# #endif(VEGAS_FOUND)

# ## ==============================================================================
# ##  PKG-CONFIG
# ## ==============================================================================

# ##Find PkgConfig which is required for VEGAS to work
# find_package(PkgConfig REQUIRED)
# if(PKG_CONFIG_FOUND)
#   message(STATUS "PkgConfig found.")
# else(NOT PKG_CONFIG_FOUND)
#   message(FATAL_ERROR "PkgConfig not found")
# endif(PKG_CONFIG_FOUND)


# ## ==============================================================================
# ##  VERITAS BANK FORMAT
# ## ==============================================================================

# pkg_search_module(VBF 0.3.4 REQUIRED VBF)
# if(VBF_FOUND)
# 	message("VBF include dirs: ${VBF_INCLUDE_DIRS}")
# 	message("VBF lib dirs: ${VBF_LIBRARIES}")
# endif(VBF_FOUND)

# ## ==============================================================================
# ##  BOOST C++ LIBRARIES
# ## ==============================================================================

# find_package(Boost)
# if(Boost_FOUND)
#   message(STATUS "BOOST ${Boost_VERSION} found!")
#   include_directories(${Boost_INCLUDE_DIR})
# else(NOT Boost_FOUND)
#   message(STATUS "BOOST not found. You need BOOST for this program to work!")
# endif(Boost_FOUND)

# ## ==============================================================================
# ##
# ##  Set directories to include in project
# ##
# ## ==============================================================================

# # file(GLOB_RECURSE vegas_headers ${VEGAS_INCLUDE_DIRS}
# #   "*.h"
# #   "*.hpp"
# #   )

# # set(SSTGATE_INCLUDE_DIRS "")
# # foreach(_headerFile ${vegas_headers})
# #     #get_filename_component(_dir ${_headerFile} PATH)
# #     get_filename_component(_dir ${_headerFile} DIRECTORY)
# #     list(APPEND SSTGATE_INCLUDE_DIRS ${_dir})
# # endforeach()
# # list(REMOVE_DUPLICATES SSTGATE_INCLUDE_DIRS)


# # foreach(_directory ${SSTGATE_INCLUDE_DIRS})
# #   message(STATUS ${_directory})
# # endforeach()

# add_definitions(-D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_OAWG -DNOHDF5 -Dnothreads)#-D_GLIBCXX_USE_CXX11_ABI=0  -D_OAWG_DIR="${CMAKE_CURRENT_BINARY_DIR}"

# if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#   add_definitions(-D__APPLE__)
# endif()

set(robast_path "$ENV{ROBAST}")
message(STATUS "ROBAST found: ${robast_path}")
include_directories("${robast_path}/include")
find_library(ROBAST_LIBRARY ROBAST HINTS ${robast_path})

include_directories(include)


message(STATUS "list dir: ${CMAKE_CURRENT_LIST_DIR}")

add_subdirectory(include)
add_subdirectory(src)

## For Doxygen documentation
#add_subdirectory(doc)

#install the version information file
install(FILES 
  ${version_outfile}
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  )
