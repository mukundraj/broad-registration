cmake_minimum_required(VERSION 3.22)

project(coordinate_mapper)

# Find Cairo
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules")
find_package(Cairo REQUIRED)

find_package(ITK REQUIRED)

include(${ITK_USE_FILE})

include_directories(./src/cpp)

add_executable(cmapper ./src/cpp/main.cpp)
# add_executable(drawer ./src/cpp/drawer.cpp)
add_executable(cmapper2 ./src/cpp/mapper2.cpp)
add_executable(img_mapper ./src/cpp/img_mapper.cpp)

target_link_libraries(cmapper ${ITK_LIBRARIES})
target_link_libraries(cmapper2 ${ITK_LIBRARIES})
target_link_libraries(img_mapper ${ITK_LIBRARIES})


# target_include_directories(drawer PRIVATE ${Cairo_INCLUDE_DIRS})
# target_link_libraries(drawer ${ITK_LIBRARIES} ${Cairo_LIBS})
