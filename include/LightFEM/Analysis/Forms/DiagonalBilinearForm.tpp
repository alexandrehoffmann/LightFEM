/*
 * DiagonalBilinearForm.tpp
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

#ifndef DIAGONAL_BILINEARFORM_TPP
#define DIAGONAL_BILINEARFORM_TPP

#include <LightFEM/Analysis/Forms/DiagonalBilinearForm.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>

template<typename Form>
DiagonalBilinearForm::DiagonalBilinearForm(const FunctionSpace *Vh, Form form) :
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	#pragma omp parallel
	{
		std::vector< double > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for collapse(2)
		for (size_t e=0;e<m_Vh->getMesh()->getNElem();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTrialFunction(e, id), m_Vh->getTestFunction(e, id));
			}
		}
		#pragma omp critical
		{
			for (size_t globId=0;globId<m_coefs.size();++globId)
			{
				m_coefs[globId] += private_coefs[globId];
			}
		}
	}
}

template<typename Form>
DiagonalBilinearForm::DiagonalBilinearForm(const FunctionSpace *Vh, Form form, MPI_Comm com) : 
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	std::vector< double > local_coefs(Vh->getNBasisFunction(), 0.0);
	
	#pragma omp parallel
	{
		MpiRange range(com, m_Vh->getMesh()->getNElem());
		
		std::vector< double > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for collapse(2)
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTrialFunction(e, id), m_Vh->getTestFunction(e, id));
			}
		}
		#pragma omp critical
		{
			for (size_t globId=0;globId<m_coefs.size();++globId)
			{
				local_coefs[globId] += private_coefs[globId];
			}
		}
	}
	MPI_Allreduce(local_coefs.data(), m_coefs.data(), m_coefs.size(), MPI_DOUBLE, MPI_SUM, com);
}

////////////////////////////////////////////////////////////////////////

template<typename Form>
CpxDiagonalBilinearForm::CpxDiagonalBilinearForm(const FunctionSpace *Vh, Form form) :
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	#pragma omp parallel
	{
		std::vector< std::complex< double > > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for collapse(2)
		for (size_t e=0;e<m_Vh->getMesh()->getNElem();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTrialFunction(e, id), m_Vh->getTestFunction(e, id));
			}
		}
		#pragma omp critical
		{
			for (size_t globId=0;globId<m_coefs.size();++globId)
			{
				m_coefs[globId] += private_coefs[globId];
			}
		}
	}
}

template<typename Form>
CpxDiagonalBilinearForm::CpxDiagonalBilinearForm(const FunctionSpace *Vh, Form form, MPI_Comm com) : 
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	std::vector< std::complex< double > > local_coefs(Vh->getNBasisFunction(), 0.0);
	
	#pragma omp parallel
	{
		MpiRange range(com, m_Vh->getMesh()->getNElem());
		
		std::vector< std::complex< double > > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for collapse(2)
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTrialFunction(e, id), m_Vh->getTestFunction(e, id));
			}
		}
		#pragma omp critical
		{
			for (size_t globId=0;globId<m_coefs.size();++globId)
			{
				local_coefs[globId] += private_coefs[globId];
			}
		}
	}
	MPI_Allreduce(local_coefs.data(), m_coefs.data(), m_coefs.size(), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, com);
}

#endif // DIAGONAL_BILINEARFORM_TPP
