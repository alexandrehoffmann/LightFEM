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

#ifndef ELEMENT_WISE_FUNCTION_OPERATORS_TPP
#define ELEMENT_WISE_FUNCTION_OPERATORS_TPP

#include <LightFEM/Expression/ElementWise/operators.hpp>

template<typename Expr>
double min(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double ret = expr[0];
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			ret = std::min(expr[Element::index2d(i,j)].eval(), ret);
		}
	}
	return ret;
}

template<typename Expr>
std::complex< double > min(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > ret = expr[0];
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			ret = std::min(expr[Element::index2d(i,j)].eval(), ret);
		}
	}
	return ret;
}

template<typename Expr>
double max(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double ret = expr[0];
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			ret = std::max(expr[Element::index2d(i,j)].eval(), ret);
		}
	}
	return ret;
}

template<typename Expr>
std::complex< double > max(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > ret = expr[0];
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			ret = std::max(expr[Element::index2d(i,j)].eval(), ret);
		}
	}
	return ret;
}


#endif // ELEMENT_WISE_FUNCTION_OPERATORS_TPP
