/*
 * VectorTrialFunction.cpp
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

#include <LightFEM/Analysis/VectorTrialFunction.hpp>

#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>

VectorTrialFunction::VectorTrialFunction(const VectorFunctionSpace* functionSpace, const Element* element, const int locId) :
	m_functionSpace(functionSpace),
	m_mesh(functionSpace->getMesh()),
	m_elem(element),
	m_locId(locId),
	m_discBasisFunction(functionSpace->getDiscBasisFunction(locId))
{
}

VectorTrialFunction::VectorTrialFunction(const VectorTrialFunction& other) :
	m_functionSpace(other.getFunctionSpace()),
	m_mesh(other.getMesh()),
	m_elem(other.getElement()),
	m_locId(other.getLocId()),
	m_discBasisFunction(other.m_discBasisFunction)
{
}

////////////////////////////////////////////////////////////////////////

GradVectorTrialFunction::GradVectorTrialFunction(const VectorTrialFunction& testFunction) :
	m_functionSpace(testFunction.getFunctionSpace()),
	m_mesh(testFunction.getMesh()),
	m_elem(testFunction.getElement()),
	m_locId(testFunction.getLocId()),
	m_discGradXiBasisFunctions(testFunction.getFunctionSpace()->getDiscGradXiBasisFunction(testFunction.getLocId()))
{
}

