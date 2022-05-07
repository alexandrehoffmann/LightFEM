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
 
#ifndef FUNCTION_OPERATORS_TPP
#define FUNCTION_OPERATORS_TPP

#include <LightFEM/Expression/Function/operators.hpp>

template<typename Expr>
double min(const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double ret = min(expr[0]);
	
	for (size_t e=1;e<expr.getMesh()->getNElem();++e)
	{
		ret = std::min(ret, min(expr[e]));
	}
	return ret;
}

template<typename Expr>
std::complex< double > min(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > ret = min(expr[0]);
	
	for (size_t e=1;e<expr.getMesh()->getNElem();++e)
	{
		ret = std::min(ret, min(expr[e]));
	}
	return ret;
}


template<typename Expr>
double max(const FunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	double ret = max(expr[0]);
	
	for (size_t e=1;e<expr.getMesh()->getNElem();++e)
	{
		ret = std::max(ret, max(expr[e]));
	}
	return ret;
}

template<typename Expr>
std::complex< double > max(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr)
{
	std::complex< double > ret = max(expr[0]);
	
	for (size_t e=1;e<expr.getMesh()->getNElem();++e)
	{
		ret = std::max(ret, max(expr[e]));
	}
	return ret;
}

#endif // FUNCTION_OPERATORS_TPP
