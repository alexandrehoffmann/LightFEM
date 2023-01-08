/*
 * ElementWiseFunction.tpp
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

#ifndef ELEMENT_WISE_FUNCTION_TPP
#define ELEMENT_WISE_FUNCTION_TPP

#include <LightFEM/Expression/ElementWise/ElementWiseFunction.hpp>

template<ExprType Type>	template<typename Expr> 
ElementWiseFunction<Type>::ElementWiseFunction(const ElementWiseFunctionExpression< Type, Expr >& expr) :
	m_values(Element::getNxiNd()), 
	m_containsTrial(expr.containsTrial()), 
	m_containsTest(expr.containsTest()), 
	m_mesh(expr.getMesh()),
	m_element(expr.getElement())
{
	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
}

template<ExprType Type>	template<typename Expr> 
ElementWiseFunction<Type>& ElementWiseFunction<Type>::operator=  (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr.getMesh();
	m_element = expr.getElement();
	
	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
ElementWiseFunction<Type>& ElementWiseFunction<Type>::operator+= (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }
	
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] += expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr> 
ElementWiseFunction<Type>& ElementWiseFunction<Type>::operator-= (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }
	
	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] -= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
ElementWiseFunction<Type>& ElementWiseFunction<Type>::operator*= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] *= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
ElementWiseFunction<Type>& ElementWiseFunction<Type>::operator/= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] /= expr[k];
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>::CpxElementWiseFunction(const ElementWiseFunctionExpression< Type, Expr >& expr) :
	m_values(Element::getNxiNd()),
	m_containsTrial(expr.containsTrial()),
	m_containsTest(expr.containsTest()),
	m_mesh(expr.getMesh()),
	m_element(expr.getElement())
{
	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>::CpxElementWiseFunction(const CpxElementWiseFunctionExpression< Type, Expr >& expr) :
	m_values(Element::getNxiNd()),
	m_containsTrial(expr.containsTrial()),
	m_containsTest(expr.containsTest()),
	m_mesh(expr.getMesh()),
	m_element(expr.getElement())
{
	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator=  (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr.getMesh();
	m_element = expr.getElement();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator=  (const CpxElementWiseFunctionExpression< Type, Expr >& expr)
{
	m_containsTrial = expr.containsTrial();
	m_containsTest = expr.containsTest();
	m_mesh = expr.getMesh();
	m_element = expr.getElement();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] = expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator+= (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] += expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator+= (const CpxElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] += expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator-= (const ElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] -= expr[k];
	}
	return *this;
}


template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator-= (const CpxElementWiseFunctionExpression< Type, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] -= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator*= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] *= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator*= (const CpxElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] *= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator/= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] /= expr[k];
	}
	return *this;
}

template<ExprType Type>	template<typename Expr>
CpxElementWiseFunction<Type>& CpxElementWiseFunction<Type>::operator/= (const CpxElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr)
{
	if (   m_mesh != expr.getMesh())    { throw std::invalid_argument("Both functions must be defined on the same mesh.");    }
	if (m_element != expr.getElement()) { throw std::invalid_argument("Both functions must be defined on the same element."); }

	m_containsTrial = m_containsTrial and expr.containsTrial();
	m_containsTest = m_containsTest and expr.containsTest();

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_values[k] /= expr[k];
	}
	return *this;
}

#endif // ELEMENT_WISE_FUNCTION_TPP
