/*
 * DiagonalBilinearForm.hpp
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

#ifndef DIAGONAL_BILINEARFORM_HPP
#define DIAGONAL_BILINEARFORM_HPP

#include <vector>
#include <complex>
#include <functional>

#include <LightFEM/Analysis/TestFunction.hpp>
#include <LightFEM/Analysis/TrialFunction.hpp>

#include <LightFEM/Expression/Function/FunctionExpression.hpp>

#include <LightFEM/Tools/MpiRange.hpp>

class FunctionSpace;

class DiagonalBilinearForm
{
public:
	template<typename Form> DiagonalBilinearForm(const FunctionSpace* Vh, Form form);
	template<typename Form> DiagonalBilinearForm(const FunctionSpace* Vh, Form form, MPI_Comm com);
	
	void setZeroOnBoundary();
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames);
	
	inline double  operator[] (const size_t i) const { return m_coefs[i]; }

	inline const double* getCoefData() const { return m_coefs.data(); }

	inline size_t getNCoefs() const { return m_coefs.size(); }
	
	inline const FunctionSpace* getTrialFunctionSpace () const { return m_Vh; } 
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; } 
private:
	const FunctionSpace*  m_Vh;
	std::vector< double > m_coefs;
};

////////////////////////////////////////////////////////////////////////

class CpxDiagonalBilinearForm
{
public:
	template<typename Form> CpxDiagonalBilinearForm(const FunctionSpace* Vh, Form form);
	template<typename Form> CpxDiagonalBilinearForm(const FunctionSpace* Vh, Form form, MPI_Comm com);
	
	void setZeroOnBoundary();
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames);
	
	inline std::complex<double> operator[] (const size_t i) const { return m_coefs[i]; }

	inline const std::complex<double>* getCoefData() const { return m_coefs.data(); }

	inline size_t getNCoefs() const { return m_coefs.size(); }
	
	inline const FunctionSpace* getTrialFunctionSpace () const { return m_Vh; } 
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; } 
private:
	const FunctionSpace*                m_Vh;
	std::vector< std::complex<double> > m_coefs;
};

#include <LightFEM/Analysis/Forms/DiagonalBilinearForm.tpp>

#endif // DIAGONAL_BILINEARFORM_HPP
