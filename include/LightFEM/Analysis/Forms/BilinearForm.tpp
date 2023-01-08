/*
 * BilinearForm.tpp
 * 
 * Copyright 2023 Alexandre Hoffmann <alexandre.hoffmann.etu@gmail.com>
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

#ifndef BILINEAR_FORM_TPP
#define BILINEAR_FORM_TPP

#include <LightFEM/Analysis/Forms/BilinearForm.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>

#include <LightFEM/Mesh/Mesh.hpp>

#include <LightFEM/LightFEM_MPI.hpp>

template<typename Form>
BilinearForm::BilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, Form form) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	std::vector< MatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for collapse(3)
	for (size_t e=0;e<m_Uh->getMesh()->getNElem();++e)
	{
		for (size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}
	initEntries(entries);
}

template<typename Form>
BilinearForm::BilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, Form form, MPI_Comm com) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }
	
	MpiRange range(com, m_Uh->getMesh()->getNElem());

	std::vector< MatrixEntry > local_entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for collapse(3)
	for (size_t e=range.begin();e<range.end();++e)
	{
		for (size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}	
	std::vector< MatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
	
	MPI_Allreduce(local_entries.data(), entries.data(), local_entries.size(), LightFEM_MPI::MPI_MATRIX_ENTRY, LightFEM_MPI::MPI_MATRIX_ENTRY_SUM, com);
	
	initEntries(entries);
}

////////////////////////////////////////////////////////////////////////

template<typename Form>
CpxBilinearForm::CpxBilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, Form form) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	std::vector< CpxMatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for collapse(3)
	for (size_t e=0;e<m_Uh->getMesh()->getNElem();++e)
	{
		for (size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}
	initEntries(entries);
}

template<typename Form>
CpxBilinearForm::CpxBilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, Form form, MPI_Comm com) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	MpiRange range(com, m_Uh->getMesh()->getNElem());

	std::vector< CpxMatrixEntry > local_entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for collapse(3)
	for (size_t e=range.begin();e<range.end();++e)
	{
		for (size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				local_entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}
	std::vector< CpxMatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
	
	MPI_Allreduce(local_entries.data(), entries.data(), local_entries.size(), LightFEM_MPI::MPI_CPX_MATRIX_ENTRY, LightFEM_MPI::MPI_CPX_MATRIX_ENTRY_SUM, com);
	
	initEntries(entries);
}

#endif // BILINEAR_FORM_TPP
