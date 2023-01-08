cmake_minimum_required(VERSION 3.0)

project(LightFEM VERSION 1.0.0 LANGUAGES CXX DESCRIPTION "Lightweigh HP finite-element library")

set(CMAKE_VERBOSE_MAKEFILE 0)

########################################################################
####                             Build                              ####
########################################################################

find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(MPI REQUIRED IMPORTED_TARGET mpi-cxx)

add_definitions(-DOMPI_SKIP_MPICXX)

set(SRCS
	src/LightFEM/LightFEM_MPI.cpp
	#
	src/LightFEM/Tools/Linesearch.cpp
	src/LightFEM/Tools/MpiRange.cpp
	#
	src/LightFEM/Mesh/Element.cpp
	src/LightFEM/Mesh/Mesh.cpp
	src/LightFEM/Mesh/GmshMesh.cpp
	src/LightFEM/Mesh/RectangularMesh.cpp
	src/LightFEM/Mesh/operators.cpp
	#
	src/LightFEM/Expression/LinAlg/RankFourTensor.cpp
	src/LightFEM/Expression/LinAlg/RankTwoTensor.cpp
	src/LightFEM/Expression/LinAlg/Matrix.cpp
	#
	src/LightFEM/Analysis/TrialFunction.cpp
	src/LightFEM/Analysis/TestFunction.cpp
	src/LightFEM/Analysis/VectorTrialFunction.cpp
	src/LightFEM/Analysis/VectorTestFunction.cpp
	#
	src/LightFEM/Analysis/Measure/ElementWiseInterpolator.cpp
	src/LightFEM/Analysis/Measure/GaussLobattoQuadrature.cpp
	src/LightFEM/Analysis/Measure/SamplingFunction.cpp
	#
	src/LightFEM/Analysis/FunctionSpace/FunctionSpace.cpp
	src/LightFEM/Analysis/FunctionSpace/PnFunctionSpace.cpp
	src/LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.cpp
	#
	src/LightFEM/Analysis/Forms/LinearForm.cpp
	src/LightFEM/Analysis/Forms/VectorLinearForm.cpp
	src/LightFEM/Analysis/Forms/BilinearForm.cpp
	src/LightFEM/Analysis/Forms/DiagonalBilinearForm.cpp
	src/LightFEM/Analysis/Forms/VectorBilinearForm.cpp
)

set(HEADERS
	include/LightFEM/LightFEM_MPI.hpp
	#
	include/LightFEM/Tools/FixedSizeLuSolver.hpp
	include/LightFEM/Tools/FixedSizeLuSolver.tpp
	include/LightFEM/Tools/FixedSizeMatrix.hpp
	include/LightFEM/Tools/FixedSizeVector.hpp
	include/LightFEM/Tools/Linesearch.hpp
	include/LightFEM/Tools/MpiRange.hpp
	#
	include/LightFEM/Mesh.hpp
	include/LightFEM/Mesh/Element.hpp
	include/LightFEM/Mesh/BoundaryElement.hpp
	include/LightFEM/Mesh/Mesh.hpp
	include/LightFEM/Mesh/GmshMesh.hpp
	include/LightFEM/Mesh/RectangularMesh.hpp
	include/LightFEM/Mesh/RectangularMesh.tpp
	include/LightFEM/Mesh/Node.hpp
	include/LightFEM/Mesh/operators.hpp
	include/LightFEM/Mesh/operators.tpp
	#
	include/LightFEM/Expression.hpp
	include/LightFEM/Expression/Traits.hpp
	#
	include/LightFEM/Expression/LinAlg/BinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/ExprType.hpp
	include/LightFEM/Expression/LinAlg/Inner.hpp
	include/LightFEM/Expression/LinAlg/Inner.tpp
	include/LightFEM/Expression/LinAlg/Norm.hpp
	include/LightFEM/Expression/LinAlg/Norm.tpp
	include/LightFEM/Expression/LinAlg/MatrixBinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/MatrixBinaryExpression.tpp
	include/LightFEM/Expression/LinAlg/MatrixExpression.hpp
	include/LightFEM/Expression/LinAlg/Matrix.hpp
	include/LightFEM/Expression/LinAlg/Matrix.tpp
	include/LightFEM/Expression/LinAlg/MatrixUnaryExpression.hpp
	include/LightFEM/Expression/LinAlg/operators.hpp
	include/LightFEM/Expression/LinAlg/Outer.hpp
	include/LightFEM/Expression/LinAlg/RankFourTensorBinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/RankFourTensorExpression.hpp
	include/LightFEM/Expression/LinAlg/RankFourTensor.hpp
	include/LightFEM/Expression/LinAlg/RankFourTensor.tpp
	include/LightFEM/Expression/LinAlg/RankFourTensorUnaryExpression.hpp
	include/LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.tpp
	include/LightFEM/Expression/LinAlg/RankTwoTensorExpression.hpp
	include/LightFEM/Expression/LinAlg/RankTwoTensor.hpp
	include/LightFEM/Expression/LinAlg/RankTwoTensor.tpp
	include/LightFEM/Expression/LinAlg/RankTwoTensorUnaryExpression.hpp
	include/LightFEM/Expression/LinAlg/RefTypeSelector.hpp
	include/LightFEM/Expression/LinAlg/ScalarBinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/ScalarExpression.hpp
	include/LightFEM/Expression/LinAlg/Scalar.hpp
	include/LightFEM/Expression/LinAlg/ScalarUnaryExpression.hpp
	include/LightFEM/Expression/LinAlg/UnaryExpression.hpp
	include/LightFEM/Expression/LinAlg/VectorBinaryExpression.hpp
	include/LightFEM/Expression/LinAlg/VectorExpression.hpp
	include/LightFEM/Expression/LinAlg/Vector.hpp
	include/LightFEM/Expression/LinAlg/VectorUnaryExpression.hpp
	#
	include/LightFEM/Expression/ElementWise/ElementWiseFunctionBinaryExpression.hpp
	include/LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp
	include/LightFEM/Expression/ElementWise/ElementWiseFunction.hpp
	include/LightFEM/Expression/ElementWise/ElementWiseFunction.tpp
	include/LightFEM/Expression/ElementWise/ElementWiseFunctionUnaryExpression.hpp
	include/LightFEM/Expression/ElementWise/ElementWiseConst.hpp
	include/LightFEM/Expression/ElementWise/ElementWiseDiff.hpp
	include/LightFEM/Expression/ElementWise/Traits.hpp
	include/LightFEM/Expression/ElementWise/operators.hpp
	include/LightFEM/Expression/ElementWise/operators.tpp
	#
	include/LightFEM/Expression/Function/FunctionBinaryExpression.hpp
	include/LightFEM/Expression/Function/FunctionExpression.hpp
	include/LightFEM/Expression/Function/Function.hpp
	include/LightFEM/Expression/Function/Function.tpp
	include/LightFEM/Expression/Function/FunctionUnaryExpression.hpp
	include/LightFEM/Expression/Function/Const.hpp
	include/LightFEM/Expression/Function/Traits.hpp
	include/LightFEM/Expression/Function/operators.hpp
	include/LightFEM/Expression/Function/operators.tpp
	#
	include/LightFEM/Expression/FiniteElementFunction/CoefBinaryOp.hpp
	include/LightFEM/Expression/FiniteElementFunction/CoefUnaryOP.hpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionBinaryExpression.hpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.hpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.tpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionUnaryExpression.hpp
	include/LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionDiffExpression.hpp
	include/LightFEM/Expression/FiniteElementFunction/Traits.hpp
	include/LightFEM/Expression/FiniteElementFunction/operators.hpp
	#
	include/LightFEM/Analysis.hpp
	include/LightFEM/Analysis/TrialFunction.hpp
	include/LightFEM/Analysis/TestFunction.hpp
	include/LightFEM/Analysis/VectorTrialFunction.hpp
	include/LightFEM/Analysis/VectorTestFunction.hpp
	include/LightFEM/Analysis/operators.hpp
	include/LightFEM/Analysis/operators.tpp
	#
	include/LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp
	include/LightFEM/Analysis/Measure/ElementWiseInterpolator.tpp
	include/LightFEM/Analysis/Measure/Measure.hpp
	include/LightFEM/Analysis/Measure/GaussLobattoQuadrature.hpp
	include/LightFEM/Analysis/Measure/SamplingFunction.hpp
	#
	include/LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp
	include/LightFEM/Analysis/FunctionSpace/PnFunctionSpace.hpp
	include/LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp
	#
	include/LightFEM/Analysis/Forms/LinearForm.hpp
	include/LightFEM/Analysis/Forms/LinearForm.tpp
	include/LightFEM/Analysis/Forms/VectorLinearForm.hpp
	include/LightFEM/Analysis/Forms/VectorLinearForm.tpp
	include/LightFEM/Analysis/Forms/BilinearForm.hpp
	include/LightFEM/Analysis/Forms/BilinearForm.tpp
	include/LightFEM/Analysis/Forms/DiagonalBilinearForm.hpp
	include/LightFEM/Analysis/Forms/DiagonalBilinearForm.tpp
	include/LightFEM/Analysis/Forms/VectorBilinearForm.hpp
	include/LightFEM/Analysis/Forms/VectorBilinearForm.tpp
)

# Declare the library target.
add_library(${PROJECT_NAME} SHARED ${SRCS} ${HEADERS})

# Configure compiler options.
target_compile_options(${PROJECT_NAME} PRIVATE -Wundef -Wvarargs -Wextra -Wall -Winit-self -Wold-style-cast -Woverloaded-virtual -Wuninitialized -Wmissing-declarations -Winit-self -flax-vector-conversions -ansi -pedantic -pedantic-errors)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

# Link dependencies.
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX PkgConfig::MPI)

# Configure the directories to search for header files.
target_include_directories(${PROJECT_NAME} PRIVATE include)

# Set the version property.
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

# Set the shared object version property to the project's major version.
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})

########################################################################
####                          Installation                          ####
########################################################################

# For access to standard installation directory variables (CMAKE_INSTALL_xDIR).
include(GNUInstallDirs)

# Set library shared object file to install.
install(TARGETS ${PROJECT_NAME}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Set library API header file to install.
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

########################################################################
####                    pkg-config Installation                     ####
########################################################################

# Create the pkg-config file from the template.
configure_file(${PROJECT_NAME}.pc.in ${PROJECT_NAME}.pc @ONLY)

# Set pkg-config file to install.
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
