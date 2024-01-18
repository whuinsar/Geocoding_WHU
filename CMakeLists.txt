cmake_minimum_required(VERSION 3.21)

project(DogGeocoding CXX)

find_package(fmt CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem serialization)
find_package(Eigen3 CONFIG REQUIRED)
find_package(TIFF REQUIRED)
find_package(gdal REQUIRED)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(geocoding)
  