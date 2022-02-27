cmake_minimum_required(VERSION 3.13)
project(PARALLEL_SKETCHES LANGUAGES C CXX)

# require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

# MPI dependencies
find_package(MPI REQUIRED)

add_executable(pkll pkll.cpp)
add_executable(preq preq.cpp)
add_executable(pudds pudds.cpp uddsketch/uddsketch.cpp)

if(TARGET MPI::MPI_CXX) # Use target if available (CMake >= 3.9)
  target_link_libraries(pkll PUBLIC MPI::MPI_CXX)
  target_link_libraries(preq PUBLIC MPI::MPI_CXX)
  target_link_libraries(pudds PUBLIC MPI::MPI_CXX)
else()
  target_compile_options(pkll PUBLIC "${MPI_CXX_COMPILE_FLAGS}")
  target_link_libraries(pkll PUBLIC "${MPI_CXX_LIBRARIES}" "${MPI_CXX_LINKFLAGS}")
  target_include_directories(pkll PUBLIC "${MPI_CXX_INCLUDE_PATH}")
  target_compile_options(preq PUBLIC "${MPI_CXX_COMPILE_FLAGS}")
  target_link_libraries(preq PUBLIC "${MPI_CXX_LIBRARIES}" "${MPI_CXX_LINKFLAGS}")
  target_include_directories(preq PUBLIC "${MPI_CXX_INCLUDE_PATH}")
  target_compile_options(pudds PUBLIC "${MPI_CXX_COMPILE_FLAGS}")
  target_link_libraries(pudds PUBLIC "${MPI_CXX_LIBRARIES}" "${MPI_CXX_LINKFLAGS}")
  target_include_directories(pudds PUBLIC "${MPI_CXX_INCLUDE_PATH}")
endif()

target_include_directories( pkll PUBLIC ${CMAKE_SOURCE_DIR}/DataSketches)
target_include_directories( pkll PUBLIC /usr/local/include)
target_include_directories( preq PUBLIC ${CMAKE_SOURCE_DIR}/DataSketches)
target_include_directories( preq PUBLIC /usr/local/include)
target_include_directories( pudds PUBLIC ${CMAKE_SOURCE_DIR}/uddsketch)
target_include_directories( pudds PUBLIC /usr/local/include)

