cmake_minimum_required(VERSION 3.20)
project(mom VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)

# find_package(GTest REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(implot CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

# add_subdirectory(dependencies)
# add_subdirectory(src)
add_executable(
  ${PROJECT_NAME}
  "./src/main.cpp"
  "./src/gui.cpp"
)
# target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)
target_link_libraries(${PROJECT_NAME}
  PRIVATE
  spdlog::spdlog
  glfw
  glad::glad
  imgui::imgui
  implot::implot
  ${OPENGL_LIBRARIES}
)