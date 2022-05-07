cmake_minimum_required(VERSION 3.0)

project(tests)

set(CMAKE_VERBOSE_MAKEFILE 0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)

# Create imported target with name PkgConfig::CONVERT
pkg_check_modules(LIGHTFEM REQUIRED IMPORTED_TARGET LightFEM)
pkg_check_modules(EIGEN3 REQUIRED IMPORTED_TARGET eigen3)
pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)

add_definitions(-DOMPI_SKIP_MPICXX)

include_directories(.)

add_executable(eigen_marine eigen_marine.cpp)
add_executable(eigen_ABC eigen_ABC.cpp)
add_executable(eigen_ABC_mass_lumping eigen_ABC_mass_lumping.cpp)
add_executable(custom_ABC_mass_lumping custom_ABC_mass_lumping.cpp)
add_executable(eigen_reverse_time eigen_reverse_time.cpp)
add_executable(eigen_wave_equation eigen_wave_equation.cpp)
add_executable(eigen_wave_equation_mass_lumping eigen_wave_equation_mass_lumping.cpp)
add_executable(cblas_wave_equation_mass_lumping cblas_wave_equation_mass_lumping.cpp)
add_executable(custom_wave_equation_mass_lumping custom_wave_equation_mass_lumping.cpp)

target_compile_options(eigen_marine PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(eigen_ABC PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(eigen_ABC_mass_lumping PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(eigen_wave_equation PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(eigen_reverse_time PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(eigen_wave_equation_mass_lumping PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(cblas_wave_equation_mass_lumping PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)
target_compile_options(custom_wave_equation_mass_lumping PRIVATE -Wundef -Wvarargs -fno-elide-constructors -pedantic-errors -ansi -Wextra -Wall -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions)

target_compile_features(eigen_marine PUBLIC cxx_std_17)
target_compile_features(eigen_ABC PUBLIC cxx_std_17)
target_compile_features(eigen_ABC_mass_lumping PUBLIC cxx_std_17)
target_compile_features(custom_ABC_mass_lumping PUBLIC cxx_std_17)
target_compile_features(eigen_reverse_time PUBLIC cxx_std_17)
target_compile_features(eigen_wave_equation PUBLIC cxx_std_17)
target_compile_features(eigen_wave_equation_mass_lumping PUBLIC cxx_std_17)
target_compile_features(cblas_wave_equation_mass_lumping PUBLIC cxx_std_17)
target_compile_features(custom_wave_equation_mass_lumping PUBLIC cxx_std_17)

target_link_libraries(eigen_marine PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(eigen_ABC PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(eigen_ABC_mass_lumping PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(custom_ABC_mass_lumping PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(eigen_reverse_time PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(eigen_wave_equation PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(eigen_wave_equation_mass_lumping PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(cblas_wave_equation_mass_lumping PUBLIC PkgConfig::GSL PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)
target_link_libraries(custom_wave_equation_mass_lumping PUBLIC PkgConfig::LIGHTFEM OpenMP::OpenMP_CXX)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/marine.geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/marine.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circle.geo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/circle.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/make_movie.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
