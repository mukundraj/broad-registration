project(coordinate_mapper)

find_package(ITK REQUIRED)

include(${ITK_USE_FILE})

include_directories(./src/cpp)

add_executable(cmapper ./src/cpp/main.cpp)

target_link_libraries(cmapper ${ITK_LIBRARIES})
