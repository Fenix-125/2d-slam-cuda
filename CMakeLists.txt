cmake_minimum_required(VERSION 3.17)
#project(2d_slam_cuda CUDA CXX)
project(2d_slam_cuda CXX)

set(CMAKE_CUDA_STANDARD 14)

include_directories(inc)

#add_executable(fast_slam_cuda main.cu)

add_executable(fast_slam main.cpp
        #        src/particle/particle.cpp inc/particle/particle.h
        #        src/utils/ScanMatcher.cpp inc/utils/ScanMatcher.h
        #        inc/slam_integrate_header.h
        src/utils/OccupancyGrid.cpp inc/utils/OccupancyGrid.h
        )

#set_target_properties(
#        fast_slam_cuda
#        PROPERTIES
#        CUDA_SEPARABLE_COMPILATION ON)
