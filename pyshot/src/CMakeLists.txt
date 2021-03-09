FIND_PACKAGE (Eigen3 3.3 REQUIRED NO_MODULE)

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
	${PROJECT_SOURCE_DIR}/src/calc_shot.cpp
 )

TARGET_LINK_LIBRARIES (shot Eigen3::Eigen ${LZ4_LIBRARIES})