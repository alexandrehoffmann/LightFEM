/*
 * ElementWiseDiff.hpp
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

#ifndef ELEMENT_WISE_DIFF_HPP
#define ELEMENT_WISE_DIFF_HPP

#include <LightFEM/Expression/Traits.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

#include <LightFEM/Expression/LinAlg/operators.hpp>

template<ExprType Type, typename Expr>
class ElementWiseDiff : public ElementWiseFunctionExpression< Type, ElementWiseDiff<Type, Expr> > {};

template<typename Expr>
class ElementWiseDiff<ExprType::SCALAR, Expr> : public ElementWiseFunctionExpression< ExprType::SCALAR, ElementWiseDiff<ExprType::SCALAR, Expr> >
{
public:
	typedef typename Traits< ElementWiseDiff<ExprType::SCALAR, Expr> >::ReturnType ReturnType;
public:
	ElementWiseDiff(Expr&& grad, const size_t dim) : m_grad(std::forward<Expr>(grad)), m_dim(dim) {}
	
	inline ReturnType operator[] (const size_t k) const { return inner(s_extractor[m_dim], m_grad[k]); }

	inline bool containsTrial() const { return m_grad.containsTrial(); }
	inline bool containsTest()  const { return m_grad.containsTest(); }

	inline const Mesh*    getMesh()    const { return m_grad.getMesh();    }
	inline const Element* getElement() const { return m_grad.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_grad;
	size_t m_dim;
private:
	static const std::array<Vector, 2> s_extractor;
};

template<typename Expr>
class ElementWiseDiff<ExprType::VECTOR, Expr> : public ElementWiseFunctionExpression< ExprType::VECTOR, ElementWiseDiff<ExprType::VECTOR, Expr> >
{
public:
	typedef typename Traits< ElementWiseDiff<ExprType::VECTOR, Expr> >::ReturnType ReturnType;
public:
	ElementWiseDiff(Expr&& grad, const size_t dim) : m_grad(std::forward<Expr>(grad)), m_dim(dim) {}
	
	inline ReturnType operator[] (const size_t k) const { return transpose(m_grad[k])*s_extractor[m_dim]; }

	inline bool containsTrial() const { return m_grad.containsTrial(); }
	inline bool containsTest()  const { return m_grad.containsTest(); }

	inline const Mesh*    getMesh()    const { return m_grad.getMesh();    }
	inline const Element* getElement() const { return m_grad.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_grad;
	size_t m_dim;
private:
	static const std::array<Vector, 2> s_extractor;
};

template<typename Expr> const std::array<Vector, 2> ElementWiseDiff<ExprType::SCALAR, Expr>::s_extractor({ Vector({1.0, 0.0}), Vector({0.0, 1.0}) });
template<typename Expr> const std::array<Vector, 2> ElementWiseDiff<ExprType::VECTOR, Expr>::s_extractor({ Vector({1.0, 0.0}), Vector({0.0, 1.0}) });

////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr> class CpxElementWiseDiff : public CpxElementWiseFunctionExpression< Type, CpxElementWiseDiff<Type, Expr> > {};

template<typename Expr>
class CpxElementWiseDiff<ExprType::SCALAR, Expr> : public CpxElementWiseFunctionExpression< ExprType::SCALAR, CpxElementWiseDiff<ExprType::SCALAR, Expr> >
{
public:
	typedef typename Traits< ElementWiseDiff<ExprType::SCALAR, Expr> >::ReturnType ReturnType;
public:
	CpxElementWiseDiff(Expr&& grad, const size_t dim) : m_grad(std::forward<Expr>(grad)), m_dim(dim) {}
	
	inline ReturnType operator[] (const size_t k) const { return inner(s_extractor[m_dim], m_grad[k]); }

	inline bool containsTrial() const { return m_grad.containsTrial(); }
	inline bool containsTest()  const { return m_grad.containsTest(); }

	inline const Mesh*    getMesh()    const { return m_grad.getMesh();    }
	inline const Element* getElement() const { return m_grad.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_grad;
	size_t m_dim;
private:
	static const std::array<Vector, 2> s_extractor;
};

template<typename Expr>
class CpxElementWiseDiff<ExprType::VECTOR, Expr> : public CpxElementWiseFunctionExpression< ExprType::VECTOR, CpxElementWiseDiff<ExprType::VECTOR, Expr> >
{
public:
	typedef typename Traits< ElementWiseDiff<ExprType::VECTOR, Expr> >::ReturnType ReturnType;
public:
	CpxElementWiseDiff(Expr&& grad, const size_t dim) : m_grad(std::forward<Expr>(grad)), m_dim(dim) {}
	
	inline ReturnType operator[] (const size_t k) const { return transpose(m_grad[k])*s_extractor[m_dim]; }

	inline bool containsTrial() const { return m_grad.containsTrial(); }
	inline bool containsTest()  const { return m_grad.containsTest(); }

	inline const Mesh*    getMesh()    const { return m_grad.getMesh();    }
	inline const Element* getElement() const { return m_grad.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_grad;
	size_t m_dim;
private:
	static const std::array<Vector, 2> s_extractor;
};

template<typename Expr> const std::array<Vector, 2> CpxElementWiseDiff<ExprType::SCALAR, Expr>::s_extractor({ Vector({1.0, 0.0}), Vector({0.0, 1.0}) });
template<typename Expr> const std::array<Vector, 2> CpxElementWiseDiff<ExprType::VECTOR, Expr>::s_extractor({ Vector({1.0, 0.0}), Vector({0.0, 1.0}) });

#endif // ELEMENT_WISE_DIFF_HPP
