/*
 * LinearForm.hpp
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

#ifndef LINEAR_FORM_HPP
#define LINEAR_FORM_HPP

#include <vector>
#include <complex>
#include <functional>

#include <LightFEM/Analysis/TestFunction.hpp>

#include <LightFEM/Expression/Function/FunctionExpression.hpp>

#include <LightFEM/Tools/MpiRange.hpp>

class FunctionSpace;

class LinearForm
{
public:
	LinearForm(const FunctionSpace* Vh, const std::function<double(const TestFunction&)>& form);
	LinearForm(const FunctionSpace* Vh, const std::function<double(const TestFunction&)>& form, MPI_Comm com);
	
	void setZeroOnBoundary();
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames);
	
	inline double  operator[] (const size_t i) const { return m_coefs[i]; }

	inline const double* getCoefData() const { return m_coefs.data(); }

	inline size_t getNCoefs() const { return m_coefs.size(); }
	
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; }
private:
	const FunctionSpace*  m_Vh;
	std::vector< double > m_coefs;
};

////////////////////////////////////////////////////////////////////////

class CpxLinearForm
{
public:
	CpxLinearForm(const FunctionSpace* Vh, const std::function<std::complex<double>(const TestFunction&)>& form);
	CpxLinearForm(const FunctionSpace* Vh, const std::function<std::complex<double>(const TestFunction&)>& form, MPI_Comm com);
	
	void setZeroOnBoundary();
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames);
	
	inline std::complex<double> operator[] (const size_t i) const { return m_coefs[i]; }

	inline const std::complex<double>* getCoefData() const { return m_coefs.data(); }

	inline size_t getNCoefs() const { return m_coefs.size(); }
	
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; }
private:
	const FunctionSpace*                m_Vh;
	std::vector< std::complex<double> > m_coefs;
};

#endif // LINEAR_FORM_HPP
