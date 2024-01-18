add_executable(eineder eineder.cxx)
target_link_libraries(eineder PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(eineder PRIVATE ${CMAKE_SOURCE_DIR}/common)

add_executable(plain plain.cxx)
target_link_libraries(plain PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(plain PRIVATE ${CMAKE_SOURCE_DIR}/common)

add_executable(ours ours.cxx)
target_link_libraries(ours PRIVATE ${CMAKE_SOURCE_DIR}/common/Common.lib Boost::filesystem
    Boost::serialization fmt::fmt cxxopts::cxxopts Eigen3::Eigen)
target_include_directories(ours PRIVATE ${CMAKE_SOURCE_DIR}/common)

