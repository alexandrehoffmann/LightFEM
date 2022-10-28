/*
 * Function.tpp
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

#ifndef FUNCTION_TPP
#define FUNCTION_TPP

#include <LightFEM/Expression/Function/Function.hpp>

template<ExprType Type>
Function<Type>::Function(const Mesh* mesh, FunctorType f) :
	m_values(mesh->getNElem()),
	m_containsTrial(false),
	m_containsTest(false),
	m_mesh(mesh)
{
	for (size_t e=0;e<m_mesh->getNElem();++e)
	{
		m_values[e].setElement(m_mesh->getElem(e));
		for (size_t i=0;i<Element::getNxi();++i) { for (size_t j=0;j<Element::getNxi();++j)
		{
			auto [x,y] = m_values[e].getElement()->getXworld(Element::get_xi(i), Element::get_xi(j));
			m_values[e][Element::index2d(i,j)] = f(x,y);
		}}
	}
}

template<ExprType Type>	template<typename Expr> 
Function<Type>::Function(const FunctionExpression< Type, Expr >& expr) :
	m_values(expr.getMesh()->getNElem()), 
	m_containsTrial(expr.containsTrial()), 
	m_containsTest(expr.containsTest()), 
	m_mesh(expr.getMesh())
{
	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
}

template<ExprType Type>	template<typename Expr> 
Function<Type>& Function<Type>::operator=  (const FunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr.getMesh();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}

	return *this;
}

template<ExprType Type>	template<typename Expr>
Function<Type>& Function<Type>::operator+= (const FunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();
	
	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr> 
Function<Type>& Function<Type>::operator-= (const FunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
Function<Type>& Function<Type>::operator*= (const FunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
Function<Type>& Function<Type>::operator/= (const FunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////

template<ExprType Type>
CpxFunction<Type>::CpxFunction(const Mesh* mesh, CpxFunctorType f) :
	m_values(mesh->getNElem()),
	m_containsTrial(false),
	m_containsTest(false),
	m_mesh(mesh)
{
	for (size_t e=0;e<m_mesh->getNElem();++e)
	{
		m_values[e].setElement(m_mesh->getElem(e));
		for (size_t i=0;i<Element::getNxi();++i) { for (size_t j=0;j<Element::getNxi();++j)
		{
			auto [x,y] = m_values[e].getElement()->getXworld(Element::get_xi(i), Element::get_xi(j));
			m_values[e][Element::index2d(i,j)] = f(x,y);
		}}
	}
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>::CpxFunction(const FunctionExpression< Type, Expr >& expr) :
	m_values(expr.getBoundaryMesh()->getNElem()),
	m_containsTrial(expr.containsTrial()),
	m_containsTest(expr.containsTest()),
	m_mesh(expr.getMesh())
{
	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>::CpxFunction(const CpxFunctionExpression< Type, Expr >& expr) :
	m_values(expr.getBoundaryMesh()->getNElem()),
	m_containsTrial(expr.containsTrial()),
	m_containsTest(expr.containsTest()),
	m_mesh(expr.getMesh())
{
	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator=  (const FunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr(expr.getMesh());

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator=  (const CpxFunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr(expr.getMesh());

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator+= (const FunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator+= (const CpxFunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator-= (const FunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}


template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator-= (const CpxFunctionExpression< Type, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator*= (const FunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator*= (const CpxFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator/= (const FunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxFunction<Type>& CpxFunction<Type>::operator/= (const CpxFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (m_mesh != expr.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		m_values[e] = expr[e];
	}
	return *this;
}

#endif // FUNCTION_TPP
