cmake_minimum_required(VERSION 3.8.2)
project(test)

if (APPLE)
	set (CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	#Bring the headers, such as Student.h into the project
	include_directories(
		/usr/local/include/voro++
		/usr/local/include/
		/home/tolga/Documents/github/cpp_analysis_codes-master/include
	)

	##Can manually add the sources using the set command as follows:
	##set(SOURCES src/mainapp.cpp src/Student.cpp)

	#However, the file(GLOB...) allows for wildcard additions:
	file(GLOB SOURCES "/home/tolga/Documents/github/cpp_analysis_codes-master/src/*.cpp")

	# set(CMAKE_CXX_FLAGS "-Wall -L/usr/local/lib -L/usr/lib -O2 -larmadillo -lgsl -lvoro++ -lmlpack -lboost_serialization")
	# set(CMAKE_CXX_FLAGS "-L/usr/local/lib")
	# set(CMAKE_CXX_FLAGS "-Wall -O2 -L/usr/local/lib -larmadillo -lgsl -lvoro++ -lmlpack -lboost_serialization")
	set(CMAKE_CXX_FLAGS "-Wall -O2 -L/usr/local/lib -larmadillo -lmlpack -lboost_serialization")

	add_executable(read_dump test_read_dump.cpp ${SOURCES})

	# target_link_libraries(read_dump voro++)
	target_link_libraries(read_dump armadillo)
	# target_link_libraries(read_dump gsl)
endif ()

if (UNIX AND NOT APPLE)
	set (CMAKE_CXX_STANDARD 11)

    # find_package(mlpack)

	#Bring the headers, such as Student.h into the project
	include_directories(
		/usr/include/gsl
		/usr/include/mlpack
		/usr/local/include/voro++
		/home/tolga/Documents/github/cpp_analysis_codes-master/include
	)

	##Can manually add the sources using the set command as follows:
	##set(SOURCES src/mainapp.cpp src/Student.cpp)

	#However, the file(GLOB...) allows for wildcard additions:
	file(GLOB SOURCES "/home/tolga/Documents/github/cpp_analysis_codes-master/src/*.cpp")

	# set(CMAKE_CXX_FLAGS "-Wall -L/usr/local/lib -L/usr/lib -O2 -larmadillo -lgsl -lvoro++")
	set(CMAKE_CXX_FLAGS "-Wall -O2")
	# set(CMAKE_CXX_FLAGS "-L/usr/local/lib")

	add_executable(test test.cpp ${SOURCES})
    target_link_libraries(test mlpack)
	target_link_libraries(test voro++)
	target_link_libraries(test armadillo)
	target_link_libraries(test gsl)
endif ()
