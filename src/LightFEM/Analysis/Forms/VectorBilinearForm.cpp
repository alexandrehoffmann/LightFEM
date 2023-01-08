/*
 * VectorBilinearForm.cpp
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

#include <LightFEM/Analysis/Forms/VectorBilinearForm.hpp>

#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>
#include <LightFEM/Mesh/Mesh.hpp>

#include <LightFEM/LightFEM_MPI.hpp>

void VectorBilinearForm::initEntries(std::vector< MatrixEntry >& entries)
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

void VectorBilinearForm::setIdentityOnBoundary()
{
	for (MatrixEntry& entry : m_entries)
	{
		if (m_Uh->isIdOnBoundary(entry.row) or m_Vh->isIdOnBoundary(entry.col))
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
}

void VectorBilinearForm::setZeroOnBoundary(const bool setTrialZero, const bool setTestZero)
{
	for (MatrixEntry& entry : m_entries)
	{
		if ((setTestZero and m_Uh->isIdOnBoundary(entry.row)) or (setTrialZero and m_Vh->isIdOnBoundary(entry.col)))
		{
			entry.value = 0.0;
		}
	}
}

void VectorBilinearForm::setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames)
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
}

void VectorBilinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames, const bool setTrialZero, const bool setTestZero)
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
}

void VectorBilinearForm::pruneNullEntries(const double tol)
{
	m_entries.erase(std::remove_if(std::begin(m_entries), std::end(m_entries), [tol](const MatrixEntry& entry) -> bool
	{
		if (std::fabs(entry.value) < tol) { return true; }
		return false; 
	}));
}

////////////////////////////////////////////////////////////////////////

void CpxVectorBilinearForm::initEntries(std::vector< CpxMatrixEntry >& entries)
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

void CpxVectorBilinearForm::setIdentityOnBoundary()
{
	for (CpxMatrixEntry& entry : m_entries)
	{
		if (m_Uh->isIdOnBoundary(entry.row) or m_Vh->isIdOnBoundary(entry.col))
		{
			entry.value = (entry.row == entry.col) ? 1.0 : 0.0;
		}
	}
}

void CpxVectorBilinearForm::setZeroOnBoundary(const bool setTrialZero, const bool setTestZero)
{
	for (CpxMatrixEntry& entry : m_entries)
	{
		if ((setTestZero and m_Uh->isIdOnBoundary(entry.row)) or (setTrialZero and m_Vh->isIdOnBoundary(entry.col)))
		{
			entry.value = 0.0;
		}
	}
}

void CpxVectorBilinearForm::setIdentityOnBoundary(std::initializer_list<std::string> boundaryNames)
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
}

void CpxVectorBilinearForm::setZeroOnBoundary(std::initializer_list<std::string> boundaryNames, const bool setTrialZero, const bool setTestZero)
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
}

void CpxVectorBilinearForm::pruneNullEntries(const double tol)
{
	m_entries.erase(std::remove_if(std::begin(m_entries), std::end(m_entries), [tol](const CpxMatrixEntry& entry) -> bool
	{
		if (std::fabs(entry.value) < tol) { return true; }
		return false; 
	}));
}

