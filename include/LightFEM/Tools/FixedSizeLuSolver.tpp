/*
 * FixedSizeLuSolver.tpp
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

#ifndef FIXED_SIZE_LU_SOLVER_TPP
#define FIXED_SIZE_LU_SOLVER_TPP

#include <LightFEM/Tools/FixedSizeLuSolver.hpp>

#include <cmath>

// c.f. https://stackoverflow.com/a/55131490/10532694

template<size_t Dim>
FixedSizeLuSolver<Dim>::FixedSizeLuSolver(const FixedSizeMatrix<Dim, Dim>& A)
{
	// initialization
	for (size_t i=0;i<Dim;++i)
	{
		m_piv[i] = i;
		for (size_t j=0;j<Dim;++j)
		{
			m_LU(i,j) = A(i,j);
		}
	}
	
	for (size_t k=0;k<Dim-1;++k)
	{
		// finding a pivot
		size_t max_row_index = k;
		double max_row_elem = std::fabs(m_LU(k,k));
		for (size_t i=k+1;i<Dim;++i)
		{
			if (std::fabs(m_LU(i,k)) > max_row_elem)
			{
				max_row_index = i;
				max_row_elem = std::fabs(m_LU(i,k));
			}
		}
		std::swap(m_piv[k], m_piv[max_row_index]);
		for (size_t j=0;j<Dim;++j) { std::swap(m_LU(k,j), m_LU(max_row_index,j)); }
		// LU facto
		for (size_t i=k+1;i<Dim;++i)
		{
			m_LU(i,k) /= m_LU(k,k);
			for (size_t j=k+1;j<Dim;++j) 
			{
				m_LU(i,j) -= m_LU(i,k)*m_LU(k,j);
			}
		}
	}
}

template<size_t Dim>
FixedSizeVector<Dim> FixedSizeLuSolver<Dim>::solve(const FixedSizeVector<Dim>& b)
{
	FixedSizeVector<Dim> Pb; for (size_t i=0;i<Dim;++i) { Pb[i] = b[m_piv[i]]; }
	
	return solve_bsub(solve_ufsub(Pb));
}

template<size_t Dim>
FixedSizeVector<Dim> FixedSizeLuSolver<Dim>::solve_ufsub(const FixedSizeVector<Dim>& b)
{
	// Unit row oriented forward substitution	
	// Solves for y, b = Ly
	
	FixedSizeVector<Dim> y(b);
	for (size_t i=0;i<Dim;++i)
	{
		for (size_t j=0;j<i;++j) 
		{
			y[i] -= m_LU(i,j)*y[j];
		}
	}
	
	return y;
}

template<size_t Dim>
FixedSizeVector<Dim> FixedSizeLuSolver<Dim>::solve_bsub(const FixedSizeVector<Dim>& y)
{
	// Row oriented backward substitution
	// Solves for x, y = Ux
	
	FixedSizeVector<Dim> x(y);
	for (int i=Dim-1;i>=0;--i)
	{
		for (size_t j=i+1;j<Dim;++j)
		{
			x[i] -= m_LU(i,j)*x[j];
		} 
		x[i] /= m_LU(i,i);
	}
	return x;
}

#endif // FIXED_SIZE_LU_SOLVER_TPP
