add_executable(eineder eineder.cxx)
target_link_libraries(eineder PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(eineder PRIVATE ${CMAKE_SOURCE_DIR}/common)

add_executable(plain plain.cxx)
target_link_libraries(plain PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(plain PRIVATE ${CMAKE_SOURCE_DIR}/common)

add_executable(dogs dogs.cxx)
target_link_libraries(dogs PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(dogs PRIVATE ${CMAKE_SOURCE_DIR}/common)

