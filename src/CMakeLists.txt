FIND_PACKAGE (Eigen3 3.3 REQUIRED NO_MODULE)

ADD_EXECUTABLE(shot
	${PROJECT_SOURCE_DIR}/include/shot_descriptor.h
	${PROJECT_SOURCE_DIR}/include/mesh.h
	${PROJECT_SOURCE_DIR}/src/shot_descriptor.cpp
	${PROJECT_SOURCE_DIR}/src/calc_shot.cpp
 )

TARGET_LINK_LIBRARIES (shot Eigen3::Eigen)