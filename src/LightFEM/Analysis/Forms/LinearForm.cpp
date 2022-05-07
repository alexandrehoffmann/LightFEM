/*
 * LinearForm.cpp
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

#include <LightFEM/Analysis/Forms/LinearForm.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>
#include <LightFEM/Mesh/Mesh.hpp>

LinearForm::LinearForm(const FunctionSpace *Vh, const std::function<double (const TestFunction &)> &form) :
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	#pragma omp parallel
	{
		std::vector< double > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for
		for (size_t e=0;e<m_Vh->getMesh()->getNElem();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTestFunction(e, id));
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

LinearForm::LinearForm(const FunctionSpace *Vh, const std::function<double (const TestFunction &)> &form, MPI_Comm com) : 
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	std::vector< double > local_coefs(Vh->getNBasisFunction(), 0.0);
	
	#pragma omp parallel
	{
		MpiRange range(com, m_Vh->getMesh()->getNElem());
		
		std::vector< double > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTestFunction(e, id));
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

void LinearForm::setZeroOnBoundary()
{
	for (size_t globId=0;globId<m_coefs.size();++globId)
	{
		if (m_Vh->isIdOnBoundary(globId))
		{
			m_coefs[globId] = 0.0;
		}
	}
}

void LinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames)
{
	const Mesh* mesh = m_Vh->getMesh();
	
	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}
	
	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	
	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (size_t locId=0;locId<m_Vh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
			}
		}
	}
	
	for (size_t globId=0;globId<m_coefs.size();++globId)
	{
		if (rowConstrained[globId])
		{
			m_coefs[globId] = 0.0;
		}
	}
}

////////////////////////////////////////////////////////////////////////

CpxLinearForm::CpxLinearForm(const FunctionSpace *Vh, const std::function<std::complex<double>(const TestFunction &)> &form) :
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	#pragma omp parallel
	{
		std::vector< std::complex< double > > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for
		for (size_t e=0;e<m_Vh->getMesh()->getNElem();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTestFunction(e, id));
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

CpxLinearForm::CpxLinearForm(const FunctionSpace *Vh, const std::function<std::complex<double>(const TestFunction &)> &form, MPI_Comm com) : 
	m_Vh(Vh),
	m_coefs(Vh->getNBasisFunction(), 0.0)
{
	std::vector< std::complex< double > > local_coefs(Vh->getNBasisFunction(), 0.0);
	
	#pragma omp parallel
	{
		MpiRange range(com, m_Vh->getMesh()->getNElem());
		
		std::vector< std::complex< double > > private_coefs(Vh->getNBasisFunction(), 0.0);
		#pragma omp for
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t id=0;id<m_Vh->getNBasisFunctionPerElement();++id)
			{
				private_coefs[m_Vh->getGlobalId(e, id)] += form(m_Vh->getTestFunction(e, id));
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

void CpxLinearForm::setZeroOnBoundary()
{
	for (size_t globId=0;globId<m_coefs.size();++globId)
	{
		if (m_Vh->isIdOnBoundary(globId))
		{
			m_coefs[globId] = 0.0;
		}
	}
}

void CpxLinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames)
{
	const Mesh* mesh = m_Vh->getMesh();
	
	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}
	
	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	
	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (size_t locId=0;locId<m_Vh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
			}
		}
	}
	
	for (size_t globId=0;globId<m_coefs.size();++globId)
	{
		if (rowConstrained[globId])
		{
			m_coefs[globId] = 0.0;
		}
	}
}
