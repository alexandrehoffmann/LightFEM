/*
 * FixedSizeLuSolver.hpp
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
 
#ifndef FIXED_SIZE_LU_SOLVER_HPP
#define FIXED_SIZE_LU_SOLVER_HPP

#include <LightFEM/Tools/FixedSizeMatrix.hpp>
#include <LightFEM/Tools/FixedSizeVector.hpp>

template<size_t Dim>
class FixedSizeLuSolver
{
public:
	FixedSizeLuSolver(const FixedSizeMatrix<Dim, Dim>& A);
	
	FixedSizeVector<Dim> solve(const FixedSizeVector<Dim>& b);
private:
	FixedSizeVector<Dim> solve_ufsub(const FixedSizeVector<Dim>& b);
	FixedSizeVector<Dim> solve_bsub(const FixedSizeVector<Dim>& y);
private:
	std::array< size_t, Dim> m_piv;
	FixedSizeMatrix<Dim, Dim> m_LU;
};

#include <LightFEM/Tools/FixedSizeLuSolver.tpp>

#endif // FIXED_SIZE_LU_SOLVER_HPP
