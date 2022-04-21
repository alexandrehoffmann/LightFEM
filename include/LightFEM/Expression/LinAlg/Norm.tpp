/*
 * Norm.tpp
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

#ifndef NORM_TPP
#define NORM_TPP

#include <LightFEM/Expression/LinAlg/Norm.hpp>

template<typename Expr>
UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr>::UnaryExpression(Expr&& expr)
{
	double sq_value = 0.0;
	for (size_t i=0;i<expr.getSize();++i) { const double expr_i = expr[i]; sq_value += expr_i*expr_i; }
	m_value = std::sqrt(sq_value);
}

template<typename Expr>
UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr>::UnaryExpression(Expr&& expr)
{
	double sq_value = 0.0;
	for (size_t i=0;i<expr.getNrows();++i) { for (size_t j=0;j<expr.getNcols();++j)
	{
		const double expr_ij = expr(i,j);
		sq_value += expr_ij*expr_ij;
	}}
	m_value = std::sqrt(sq_value);
}

template<typename Expr>
UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr>::UnaryExpression(Expr&& expr)
{
	const std::array<size_t, 2> shape = expr.getShape();

	double sq_value = 0.0;
	for (size_t i=0;i<shape[0];++i) { for (size_t j=0;j<shape[1];++j)
	{
		const double expr_ij = expr(i,j);
		sq_value += expr_ij*expr_ij;
	}}
	m_value = std::sqrt(sq_value);
}

////////////////////////////////////////////////////////////////////////

template<typename Expr>
CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr>::CpxUnaryExpression(Expr&& expr)
{
	std::complex< double > sq_value = 0.0;
	for (size_t i=0;i<expr.getSize();++i) { const double expr_i = expr[i]; sq_value += expr_i*std::conj(expr_i); }
	m_value = std::sqrt(sq_value).real();
}

template<typename Expr>
CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr>::CpxUnaryExpression(Expr&& expr)
{
	std::complex< double > sq_value = 0.0;
	for (size_t i=0;i<expr.getNrows();++i) { for (size_t j=0;j<expr.getNcols();++j)
	{
		const double expr_ij = expr(i,j);
		sq_value += expr_ij*std::conj(expr_ij);
	}}
	m_value = std::sqrt(sq_value).real();
}

template<typename Expr>
CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr>::CpxUnaryExpression(Expr&& expr)
{
	const std::array<size_t, 2> shape = expr.getShape();

	std::complex< double > sq_value = 0.0;
	for (size_t i=0;i<shape[0];++i) { for (size_t j=0;j<shape[1];++j)
	{
		const double expr_ij = expr(i,j);
		sq_value += expr_ij*std::conj(expr_ij);
	}}
	m_value = std::sqrt(sq_value).real();
}

#endif // NORM_TPP
