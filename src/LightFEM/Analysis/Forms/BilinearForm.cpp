/*
 * BilinearForm.cpp
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

#include <LightFEM/Analysis/Forms/BilinearForm.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>
#include <LightFEM/Mesh/Mesh.hpp>

#include <LightFEM/LightFEM_MPI.hpp>

BilinearForm::BilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, const std::function<double (const TrialFunction &, const TestFunction &)> &form) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	std::vector< MatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for
	for (std::size_t e=0;e<m_Uh->getMesh()->getNElem();++e)
	{
		for (std::size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (std::size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}
	initEntries(entries);
}

BilinearForm::BilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, const std::function<double (const TrialFunction &, const TestFunction &)> &form, MPI_Comm com) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }
	
	MpiRange range(com, m_Uh->getMesh()->getNElem());

	std::vector< MatrixEntry > local_entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for
	for (std::size_t e=range.begin();e<range.end();++e)
	{
		for (std::size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (std::size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
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
	pruneNullEntries();
}

void BilinearForm::initEntries(std::vector< MatrixEntry >& entries)
{
	std::sort(std::begin(entries), std::end(entries), [](const MatrixEntry& lhs, const MatrixEntry& rhs) -> bool
	{
		return std::tie(lhs.row, lhs.col, lhs.value) < std::tie(rhs.row, rhs.col, rhs.value);
	});

	m_entries.reserve(entries.size());

	m_entries.push_back(entries[0]);
	for (size_t i=1;i<entries.size();++i)
	{
		if (m_entries.back().row == entries[i].row and m_entries.back().col == entries[i].col)
		{
			m_entries.back().value += entries[i].value;
		}
		else
		{
			m_entries.push_back(entries[i]);
		}
	}
}

void BilinearForm::setIdentityOnBoundary()
{
	for (MatrixEntry& entry : m_entries)
	{
		if (m_Uh->isIdOnBoundary(entry.row) or m_Vh->isIdOnBoundary(entry.col))
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
	pruneNullEntries();
}

void BilinearForm::setZeroOnBoundary(const bool setTrialZero, const bool setTestZero)
{
	for (MatrixEntry& entry : m_entries)
	{
		if ((setTestZero and m_Uh->isIdOnBoundary(entry.row)) or (setTrialZero and m_Vh->isIdOnBoundary(entry.col)))
		{
			entry.value = 0.0;
		}
	}
	pruneNullEntries();
}

void BilinearForm::setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames)
{
	const Mesh* mesh = m_Uh->getMesh();

	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}

	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	std::vector< bool > colConstrained(m_Uh->getNBasisFunction(), false);

	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					colConstrained[globId] = colConstrained[globId] or m_Uh->isIdOnBoundary(globId);
				}
			}
		}
	}

	for (MatrixEntry& entry : m_entries)
	{
		if (rowConstrained[entry.row] or colConstrained[entry.col])
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
	pruneNullEntries();
}

void BilinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames, const bool setTrialZero, const bool setTestZero)
{
	const Mesh* mesh = m_Uh->getMesh();

	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}

	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	std::vector< bool > colConstrained(m_Uh->getNBasisFunction(), false);

	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					colConstrained[globId] = colConstrained[globId] or m_Uh->isIdOnBoundary(globId);
				}
			}
		}
	}

	for (MatrixEntry& entry : m_entries)
	{
		if ((setTestZero and rowConstrained[entry.row]) or (setTrialZero and colConstrained[entry.col]))
		{
			entry.value = 0.0;
		}
	}
	pruneNullEntries();
}

void BilinearForm::pruneNullEntries()
{
	m_entries.erase(std::remove_if(std::begin(m_entries), std::end(m_entries), [](const MatrixEntry& entry) -> bool
	{
		if (std::fabs(entry.value) < 1.0e-15) { return true; }
		return false; 
	}));
}

////////////////////////////////////////////////////////////////////////

CpxBilinearForm::CpxBilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, const std::function<std::complex<double>(const TrialFunction &, const TestFunction &)> &form) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	std::vector< CpxMatrixEntry > entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for
	for (std::size_t e=0;e<m_Uh->getMesh()->getNElem();++e)
	{
		for (std::size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (std::size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
			{
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].row = m_Vh->getGlobalId(e, i);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].col = m_Uh->getGlobalId(e, j);
				entries[j + m_Uh->getNBasisFunctionPerElement()*(i + e*m_Vh->getNBasisFunctionPerElement())].value = form(m_Uh->getTrialFunction(e, j), m_Vh->getTestFunction(e, i));
			}
		}
	}
	initEntries(entries);
	pruneNullEntries();
}

CpxBilinearForm::CpxBilinearForm(const FunctionSpace *Uh, const FunctionSpace *Vh, const std::function<std::complex<double>(const TrialFunction &, const TestFunction &)> &form, MPI_Comm com) :
	m_Uh(Uh),
	m_Vh(Vh)
{
	if (Uh->getMesh() != Vh->getMesh()) { throw std::invalid_argument("The two function space must be defined on the same mesh"); }

	MpiRange range(com, m_Uh->getMesh()->getNElem());

	std::vector< CpxMatrixEntry > local_entries(m_Uh->getMesh()->getNElem()*m_Vh->getNBasisFunctionPerElement()*m_Uh->getNBasisFunctionPerElement());
#pragma omp parallel for
	for (std::size_t e=range.begin();e<range.end();++e)
	{
		for (std::size_t i=0;i<m_Vh->getNBasisFunctionPerElement();++i)
		{
			for (std::size_t j=0;j<m_Uh->getNBasisFunctionPerElement();++j)
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
	pruneNullEntries();
}

void CpxBilinearForm::initEntries(std::vector< CpxMatrixEntry >& entries)
{
	std::sort(std::begin(entries), std::end(entries), [](const CpxMatrixEntry& lhs, const CpxMatrixEntry& rhs) -> bool
	{
		return std::make_tuple(lhs.row, lhs.col, lhs.value.real(), lhs.value.imag()) < std::make_tuple(rhs.row, rhs.col, rhs.value.real(), rhs.value.imag());
	});
	m_entries.reserve(entries.size());

	m_entries.push_back(entries[0]);
	for (size_t i=1;i<entries.size();++i)
	{
		if (m_entries.back().row == entries[i].row and m_entries.back().col == entries[i].col)
		{
			m_entries.back().value += entries[i].value;
		}
		else
		{
			m_entries.push_back(entries[i]);
		}
	}
}

void CpxBilinearForm::setIdentityOnBoundary()
{
	for (CpxMatrixEntry& entry : m_entries)
	{
		if (m_Uh->isIdOnBoundary(entry.row) or m_Vh->isIdOnBoundary(entry.col))
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
	pruneNullEntries();
}

void CpxBilinearForm::setZeroOnBoundary(const bool setTrialZero, const bool setTestZero)
{
	for (CpxMatrixEntry& entry : m_entries)
	{
		if ((setTestZero and m_Uh->isIdOnBoundary(entry.row)) or (setTrialZero and m_Vh->isIdOnBoundary(entry.col)))
		{
			entry.value = 0.0;
		}
	}
	pruneNullEntries();
}

void CpxBilinearForm::setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames)
{
	const Mesh* mesh = m_Uh->getMesh();

	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}

	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	std::vector< bool > colConstrained(m_Uh->getNBasisFunction(), false);

	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					colConstrained[globId] = colConstrained[globId] or m_Uh->isIdOnBoundary(globId);
				}
			}
		}
	}

	for (CpxMatrixEntry& entry : m_entries)
	{
		if (rowConstrained[entry.row] or colConstrained[entry.col])
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
	pruneNullEntries();
}

void CpxBilinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames, const bool setTrialZero, const bool setTestZero)
{
	const Mesh* mesh = m_Uh->getMesh();

	std::vector< size_t > ids;
	for (const std::string& boundaryName : boundaryNames)
	{
		ids.push_back(mesh->getBoundaryId(boundaryName));
	}

	std::vector< bool > rowConstrained(m_Vh->getNBasisFunction(), false);
	std::vector< bool > colConstrained(m_Uh->getNBasisFunction(), false);

	for (size_t be=0;be<mesh->getNBoundaryElem();++be)
	{
		for (const size_t id : ids)
		{
			if (mesh->getBoundaryElem(be)->isInDomain(id))
			{
				const size_t e = mesh->getElemIdFromBoundaryElemId(be);
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					rowConstrained[globId] = rowConstrained[globId] or m_Vh->isIdOnBoundary(globId);
				}
				for (std::size_t locId=0;locId<m_Uh->getNBasisFunctionPerElement();++locId)
				{
					const size_t globId = m_Vh->getGlobalId(e, locId);
					colConstrained[globId] = colConstrained[globId] or m_Uh->isIdOnBoundary(globId);
				}
			}
		}
	}

	for (CpxMatrixEntry& entry : m_entries)
	{
		if ((setTestZero and rowConstrained[entry.row]) or (setTrialZero and colConstrained[entry.col]))
		{
			entry.value = 0.0;
		}
	}
	pruneNullEntries();
}

void CpxBilinearForm::pruneNullEntries()
{
	m_entries.erase(std::remove_if(std::begin(m_entries), std::end(m_entries), [](const CpxMatrixEntry& entry) -> bool
	{
		if (std::fabs(entry.value) < 1.0e-15) { return true; }
		return false; 
	}));
}
