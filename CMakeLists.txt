CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(shot VERSION 0.1.0)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
SET(CMAKE_BUILD_TYPE Debug )

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
LINK_DIRECTORIES( ${PROJECT_SOURCE_DIR}/lib)

FIND_PACKAGE (Eigen3 3.3 REQUIRED NO_MODULE)
IF(Eigen3)
    MESSAGE(STATUS "Found: ${Eigen3::Eigen}")
ENDIF()
# Since flann 1.8.5, LZ4 is vendored. Building against this version fails.
# This explicitly states to use the system's version of LZ4 instead.
# See: https://github.com/mariusmuja/flann/issues/384#issuecomment-668782051
MESSAGE(STATUS "Finding LZ4 libraries")
FIND_LIBRARY(LZ4_LIBRARIES NAMES liblz4.so HINTS "lz4/lib")
IF(LZ4_LIBRARIES)
    MESSAGE(STATUS "Found: ${LZ4_LIBRARIES}")
    ELSE()
    MESSAGE(STATUS "Not found: ${LZ4_LIBRARIES}")
    MESSAGE(FATAL_ERROR "Cannot find required LZ4 libraries")
ENDIF()

ADD_EXECUTABLE(shot
	${PROJECT_SOURCE_DIR}/include/shot_descriptor.h
	${PROJECT_SOURCE_DIR}/include/mesh.h
	${PROJECT_SOURCE_DIR}/src/shot_descriptor.cpp
	${PROJECT_SOURCE_DIR}/src/main.cpp
 )

TARGET_LINK_LIBRARIES (shot Eigen3::Eigen ${LZ4_LIBRARIES})