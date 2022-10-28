/*
 * operators.tpp
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

#ifndef MESH_OPERATORS_TPP
#define MESH_OPERATORS_TPP

#include <fstream>

#include <LightFEM/Mesh/operators.hpp>
#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp>

////////////////////////////////////////////////////////////////////////
////                        print functions                         ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
void printFunction(const std::string& fname, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::ofstream out(fname);

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		const Element* element = expr.getMesh()->getElem(e);

		for (size_t i=0;i<Element::getNxi();++i)
		{
			for (size_t j=0;j<Element::getNxi();++j)
			{
				out << element->getXworld(Element::get_xi(i), Element::get_xi(j)).x << " " << element->getXworld(Element::get_xi(i), Element::get_xi(j)).y << " " << expr[e][Element::index2d(i,j)] << std::endl;
			}
			out << std::endl;
		}
		out << std::endl;
	}
}

template<typename Expr>
void printFunction(const std::string& fname, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::ofstream out(fname);

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		const Element* element = expr.getMesh()->getElem(e);

		for (size_t i=0;i<Element::getNxi();++i)
		{
			for (size_t j=0;j<Element::getNxi();++j)
			{
				out << element->getXworld(Element::get_xi(i), Element::get_xi(j)).x << " " << element->getXworld(Element::get_xi(i), Element::get_xi(j)).y << " " << expr[e][Element::index2d(i,j)].real() << " " << expr[e][Element::index2d(i,j)].imag() << std::endl;
			}
			out << std::endl;
		}
		out << std::endl;
	}
}


template<typename Expr>
void printFunctionCoarse(const std::string& fname, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::ofstream out(fname);

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		const Element* element = expr.getMesh()->getElem(e);

		for (const size_t i : {size_t(0), Element::getNxi()-1})
		{
			for (const size_t j : {size_t(0), Element::getNxi()-1})
			{
				out << element->getXworld(Element::get_xi(i), Element::get_xi(j)).x << " " << element->getXworld(Element::get_xi(i), Element::get_xi(j)).y << " " << expr[e][Element::index2d(i,j)] << std::endl;
			}
			out << std::endl;
		}
		out << std::endl;
	}
}

template<typename Expr>
void printFunctionCoarse(const std::string& fname, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::ofstream out(fname);

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		const Element* element = expr.getMesh()->getElem(e);

		for (const size_t i : {size_t(0), Element::getNxi()-1})
		{
			for (const size_t j : {size_t(0), Element::getNxi()-1})
			{
				out << element->getXworld(Element::get_xi(i), Element::get_xi(j)).x << " " << element->getXworld(Element::get_xi(i), Element::get_xi(j)).y << " " << expr[e][Element::index2d(i,j)].real() << " " << expr[e][Element::index2d(i,j)].imag() << std::endl;
			}
			out << std::endl;
		}
		out << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////
////                     integral on an element                     ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Element& element, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double value = 0.0;

	for (size_t i=1;i<Element::getNxi()-1;++i) { for (size_t j=1;j<Element::getNxi()-1;++j)
	{
		value += expr[Element::index2d(i,j)]*element.getWeight(i, j);
	}}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element& element, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > value = 0.0;

	for (size_t i=1;i<Element::getNxi()-1;++i) { for (size_t j=1;j<Element::getNxi()-1;++j)
	{
		value += expr[Element::index2d(i,j)]*element.getWeight(i, j);
	}}

	return value;
}

template<typename Expr>
double integral(const Element&, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< NodeAndWeight >& nodeAndWeights)
{
	ElementWiseInterpolator interp_expr(expr);

	double value = 0.0;

	for (const auto& [xi1, xi2, wi] : nodeAndWeights)
	{
		value += interp_expr(xi1, xi2)*wi;
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element&, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxNodeAndWeight >& nodeAndWeights)
{
	ElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	for (const auto& [xi1, xi2, wi] : nodeAndWeights)
	{
		value += interp_expr(xi1, xi2)*wi;
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element&, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< NodeAndWeight >& nodeAndWeights)
{
	CpxElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	for (const auto& [xi1, xi2, wi] : nodeAndWeights)
	{
		value += interp_expr(xi1, xi2)*wi;
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element&, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxNodeAndWeight >& nodeAndWeights)
{
	CpxElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	for (const auto& [xi1, xi2, wi] : nodeAndWeights)
	{
		value += interp_expr(xi1, xi2)*wi;
	}

	return value;
}

////////////////////////////////////////////////////////////////////////
////             integral on the boundary of an element             ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const size_t j = Element::getNxi()-1;
		for (size_t i=1;i<Element::getNxi()-1;++i)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::TOP, i);
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const size_t j = 0;
		for (size_t i=1;i<Element::getNxi()-1;++i)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::BOTTTOM, i);
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const size_t i = 0;
		for (size_t j=1;j<Element::getNxi()-1;++j)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::LEFT, j);
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const size_t i = Element::getNxi()-1;
		for (size_t j=1;j<Element::getNxi()-1;++j)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::RIGHT, j);
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const size_t j = Element::getNxi()-1;
		for (size_t i=1;i<Element::getNxi()-1;++i)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::TOP, i);
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const size_t j = 0;
		for (size_t i=1;i<Element::getNxi()-1;++i)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::BOTTTOM, i);
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const size_t i = 0;
		for (size_t j=1;j<Element::getNxi()-1;++j)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::LEFT, j);
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const size_t i = Element::getNxi()-1;
		for (size_t j=1;j<Element::getNxi()-1;++j)
		{
			value += expr[Element::index2d(i,j)]*element.getWeight(Element::Boundary::RIGHT, j);
		}
	}

	return value;
}

template<typename Expr>
double integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< BoundaryNodeAndWeight >& nodeAndWeights)
{
	ElementWiseInterpolator interp_expr(expr);

	double value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const double xi2 = Element::get_xi(Element::getNxi()-1);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.topWi;
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const double xi2 = Element::get_xi(0);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.bottomWi;
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const double xi1 = Element::get_xi(0);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.leftWi;
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const double xi1 = Element::get_xi(Element::getNxi()-1);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.rightWi;
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxBoundaryNodeAndWeight >& nodeAndWeights)
{
	ElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const double xi2 = Element::get_xi(Element::getNxi()-1);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.topWi;
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const double xi2 = Element::get_xi(0);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.bottomWi;
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const double xi1 = Element::get_xi(0);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.leftWi;
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const double xi1 = Element::get_xi(Element::getNxi()-1);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.rightWi;
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< BoundaryNodeAndWeight >& nodeAndWeights)
{
	CpxElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const double xi2 = Element::get_xi(Element::getNxi()-1);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.topWi;
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const double xi2 = Element::get_xi(0);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.bottomWi;
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const double xi1 = Element::get_xi(0);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.leftWi;
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const double xi1 = Element::get_xi(Element::getNxi()-1);
		for (const BoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.rightWi;
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > integral(const Element& element, const Element::Boundary& boundary, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const std::vector< CpxBoundaryNodeAndWeight >& nodeAndWeights)
{
	CpxElementWiseInterpolator interp_expr(expr);

	std::complex< double > value = 0.0;

	if (boundary == Element::Boundary::TOP)
	{
		const double xi2 = Element::get_xi(Element::getNxi()-1);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.topWi;
		}
	}
	else if (boundary == Element::Boundary::BOTTTOM)
	{
		const double xi2 = Element::get_xi(0);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(nodeAndWeight.t, xi2)*nodeAndWeight.bottomWi;
		}
	}
	else if (boundary == Element::Boundary::LEFT)
	{
		const double xi1 = Element::get_xi(0);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.leftWi;
		}
	}
	else if (boundary == Element::Boundary::RIGHT)
	{
		const double xi1 = Element::get_xi(Element::getNxi()-1);
		for (const CpxBoundaryNodeAndWeight& nodeAndWeight : nodeAndWeights)
		{
			value += interp_expr(xi1, nodeAndWeight.t)*nodeAndWeight.rightWi;
		}
	}

	return value;
}


////////////////////////////////////////////////////////////////////////
////                       integral on a mesh                       ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e]); }
	return value;
}

template<typename Expr>
double integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr);
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e]); }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr);
}

template<typename Expr>
double integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); }
	return value;
}

template<typename Expr>
double integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e));
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e));
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e));
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e));
}

////////////////////////////////////////////////////////////////////////
////                    integral on a subdomain                     ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) {	value += integral(*mesh.getElem(e), expr[e]); } }
	return value;
}

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) { return integral(*mesh.getElem(e), expr); }
	return 0.0;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) { value += integral(*mesh.getElem(e), expr[e]); } }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) { return integral(*mesh.getElem(e), expr); }
	return 0.0;
}

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	double value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); } }
	return value;
}

template<typename Expr>
double integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) { return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e)); }
	return 0.0;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); } }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) { return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e)); }
	return 0.0;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); } }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) { return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e)); }
	return 0.0;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	if (measure.getMesh() != &mesh or expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getDomainId(domainName);

	std::complex< double > value = 0.0;
	for (size_t e=0;e<mesh.getNElem();++e) { if (mesh.getElem(e)->isInDomain(id)) { value += integral(*mesh.getElem(e), expr[e], measure.getNodesAndWeights(e)); } }
	return value;
}

template<typename Expr>
std::complex< double > integral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& measure)
{
	const size_t e = mesh.getElemId(expr.getElement());
	if (e == mesh.getNElem()) { return 0.0; } // expr not defined on mesh

	const int id = mesh.getDomainId(domainName);

	if (mesh.getElem(e)->isInDomain(id)) {return integral(*mesh.getElem(e), expr, measure.getNodesAndWeights(e)); }
	return 0.0;
}

////////////////////////////////////////////////////////////////////////
////               integral on the boundary of a mesh               ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	double value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[mesh.getElemIdFromBoundaryElemId(be)]);
	}

	return value;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr);
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[mesh.getElemIdFromBoundaryElemId(be)]);
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr);
		}
	}
	return 0.0;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	double value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		const size_t e = mesh.getElemIdFromBoundaryElemId(be);
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
	}

	return value;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		const size_t e = mesh.getElemIdFromBoundaryElemId(be);
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		const size_t e = mesh.getElemIdFromBoundaryElemId(be);
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	std::complex< double > value = 0.0;

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		const size_t e = mesh.getElemIdFromBoundaryElemId(be);
		value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement())
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

////////////////////////////////////////////////////////////////////////
////             integral on a boundary of a subdomain              ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	double value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) {  value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[mesh.getElemIdFromBoundaryElemId(be)]); }
	}

	return value;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr);
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	std::complex< double > value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) {  value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[mesh.getElemIdFromBoundaryElemId(be)]); }
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr);
		}
	}
	return 0.0;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	double value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) 
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
		}
	}

	return value;
}

template<typename Expr>
double boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const FunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	std::complex< double > value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) 
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	std::complex< double > value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) 
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[e], sampling.getBoundaryNodesAndWeights(e));
		}
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const Measure& sampling)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			const size_t e = mesh.getElemIdFromBoundaryElemId(be);
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling.getBoundaryNodesAndWeights(e));
		}
	}
	return 0.0;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	if (expr.getMesh() != &mesh) { return 0.0; }

	const int id = mesh.getBoundaryId(domainName);

	std::complex< double > value = 0.0;
	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (mesh.getBoundaryElem(be)->isInDomain(id)) {  value += integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr[mesh.getElemIdFromBoundaryElemId(be)], sampling); }
	}

	return value;
}

template<typename Expr>
std::complex< double > boundaryIntegral(const Mesh& mesh, std::string domainName, const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr, const CpxMeasure& sampling)
{
	const int id = mesh.getBoundaryId(domainName);

	for (size_t be=0;be<mesh.getNBoundaryElem();++be)
	{
		if (expr.getElement() == mesh.getBoundaryElem(be)->getParentElement() and mesh.getBoundaryElem(be)->isInDomain(id))
		{
			return integral(*mesh.getBoundaryElem(be)->getParentElement(), mesh.getBoundaryElem(be)->getBoundary(), expr, sampling);
		}
	}
	return 0.0;
}

#endif // MESH_OPERATORS_TPP
