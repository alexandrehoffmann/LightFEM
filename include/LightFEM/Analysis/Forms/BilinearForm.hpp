/*
 * BilinearForm.hpp
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

#ifndef BILINEAR_FORM_HPP
#define BILINEAR_FORM_HPP

#include <vector>
#include <complex>
#include <functional>

#include <LightFEM/Analysis/TestFunction.hpp>
#include <LightFEM/Analysis/TrialFunction.hpp>

#include <LightFEM/Analysis/Forms/MatrixEntry.hpp>

#include <LightFEM/Tools/MpiRange.hpp>

class FunctionSpace;

class BilinearForm
{
public:
	template<typename Form> BilinearForm(const FunctionSpace* Uh, const FunctionSpace* Vh, Form form);
	template<typename Form> BilinearForm(const FunctionSpace* Uh, const FunctionSpace* Vh, Form form, MPI_Comm com);
	BilinearForm(const BilinearForm& other) : m_Uh(other.m_Uh), m_Vh(other.m_Vh), m_entries(other.m_entries) {}
	
	void setIdentityOnBoundary();
	void setZeroOnBoundary(const bool setTrialZero=true, const bool setTestZero=true);
	
	void setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames);
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames, const bool setTrialZero=true, const bool setTestZero=true);
	
	void pruneNullEntries(const double tol=1.0e-14);
	
	inline const MatrixEntry* getCoefData() const { return m_entries.data(); }
	inline const MatrixEntry& getCoef(const size_t i) const { return m_entries[i]; }
	inline size_t getNCoefs() const { return m_entries.size(); }
	
	inline const FunctionSpace* getTrialFunctionSpace () const { return m_Uh; } 
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; } 
private:
	void initEntries(std::vector< MatrixEntry >& entries);
private:
	const FunctionSpace*       m_Uh;
	const FunctionSpace*       m_Vh;
	std::vector< MatrixEntry > m_entries;
};

////////////////////////////////////////////////////////////////////////

class CpxBilinearForm
{
public:
	template<typename Form> CpxBilinearForm(const FunctionSpace* Uh, const FunctionSpace* Vh, Form form);
	template<typename Form> CpxBilinearForm(const FunctionSpace* Uh, const FunctionSpace* Vh, Form form, MPI_Comm com);
	CpxBilinearForm(const CpxBilinearForm& other) : m_Uh(other.m_Uh), m_Vh(other.m_Vh), m_entries(other.m_entries) {}
	
	void setIdentityOnBoundary();
	void setZeroOnBoundary(const bool setTrialZero=true, const bool setTestZero=true);
	
	void setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames);
	void setZeroOnBoundary(std::initializer_list< std::string > boundaryNames, const bool setTrialZero=true, const bool setTestZero=true);
	
	void pruneNullEntries(const double tol=1.0e-14);
	
	inline const CpxMatrixEntry* getCoefData() const { return m_entries.data(); }
	inline const CpxMatrixEntry& getCoef(const size_t i) const { return m_entries[i]; }
	inline size_t getNCoefs() const { return m_entries.size(); }
	
	inline const FunctionSpace* getTrialFunctionSpace () const { return m_Uh; } 
	inline const FunctionSpace* getTestFunctionSpace  () const { return m_Vh; } 
private:
	void initEntries(std::vector< CpxMatrixEntry >& entries);
private:
	const FunctionSpace*          m_Uh;
	const FunctionSpace*          m_Vh;
	std::vector< CpxMatrixEntry > m_entries;
};

#include <LightFEM/Analysis/Forms/BilinearForm.tpp>

#endif // BILINEAR_FORM_HPP
