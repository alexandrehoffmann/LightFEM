/*
 * operators.hpp
 * 
 * Copyright 2022 Alexandre Hoffmann <alexandre.hoffmann.etu@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef MESH_OPERATORS_HPP
#define MESH_OPERATORS_HPP

#include <LightFEM/Mesh/Element.hpp>
#include <LightFEM/Mesh/Mesh.hpp>

#include <LightFEM/Analysis/Measure/Measure.hpp>

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

////////////////////////////////////////////////////////////////////////
////                        print functions                         ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
void printFunction(const std::string& fname, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
void printFunction(const std::string& fname, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
void printFunctionCoarse(const std::string& fname, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
void printFunctionCoarse(const std::string& fname, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

////////////////////////////////////////////////////////////////////////
////                     integral on an element                     ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Element& element, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Element& element, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Element& element, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< NodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxNodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< NodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxNodeAndWeight >& nodeAndWeights);

////////////////////////////////////////////////////////////////////////
////             integral on the boundary of an element             ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< BoundaryNodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxBoundaryNodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< BoundaryNodeAndWeight >& nodeAndWeights);

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxBoundaryNodeAndWeight >& nodeAndWeights);

////////////////////////////////////////////////////////////////////////
////                       integral on a mesh                       ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
double integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

////////////////////////////////////////////////////////////////////////
////                    integral on a subdomain                     ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

////////////////////////////////////////////////////////////////////////
////               integral on the boundary of a mesh               ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

////////////////////////////////////////////////////////////////////////
////             integral on a boundary of a subdomain              ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling);

#include <LightFEM/Mesh/operators.tpp>

#endif // MESH_OPERATORS_HPP
