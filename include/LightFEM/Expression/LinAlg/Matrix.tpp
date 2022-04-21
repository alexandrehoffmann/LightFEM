/*
 * Matrix.tpp
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

#ifndef MATRIX_TPP
#define MATRIX_TPP

template<typename Expr>
Matrix::Matrix(const MatrixExpression< Expr >& expr) : 
	m_nrows(expr.getNrows()), 
	m_ncols(expr.getNcols()), 
	m_core(m_nrows*m_ncols) 
{ 
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{
		m_core[flatIndex(i,j)] =  expr(i,j); 
	}} 
}

template<typename Expr>
Matrix& Matrix::operator= (const MatrixExpression< Expr >& expr)
{
	m_nrows = expr.getNrows();
	m_ncols = expr.getNcols(); 
	m_core.resize(m_nrows*m_ncols);
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] =  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
Matrix& Matrix::operator+= (const MatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] +=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
Matrix& Matrix::operator-= (const MatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] -=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
Matrix& Matrix::operator*= (const ScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] *=  expr.eval(); 
	}}
	
	return *this;
}

template<typename Expr>
Matrix& Matrix::operator/= (const ScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); };
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] /=  expr.eval(); 
	}}
	
	return *this;
}

////////////////////////////////////////////////////////////////////////

template<typename Expr>
CpxMatrix& CpxMatrix::operator= (const MatrixExpression< Expr >& expr)
{
	m_nrows = expr.getNrows();
	m_ncols = expr.getNcols(); 
	m_core.resize(m_nrows*m_ncols);
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] =  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator= (const CpxMatrixExpression< Expr >& expr)
{
	m_nrows = expr.getNrows();
	m_ncols = expr.getNcols(); 
	m_core.resize(m_nrows*m_ncols);
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] =  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator+= (const MatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] +=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator+= (const CpxMatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] +=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator-= (const MatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] -=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator-= (const CpxMatrixExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] -=  expr(i,j); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator*= (const ScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] *=  expr.eval(); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator*= (const CpxScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] *=  expr.eval(); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator/= (const ScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] /=  expr.eval(); 
	}}
	
	return *this;
}

template<typename Expr>
CpxMatrix& CpxMatrix::operator/= (const CpxScalarExpression< Expr >& expr)
{
	if (getNrows() != expr.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
	if (getNcols() != expr.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); }
	
	for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) 
	{ 
		m_core[flatIndex(i,j)] /=  expr.eval(); 
	}}
	
	return *this;
}

#endif  // MATRIX_TPP
