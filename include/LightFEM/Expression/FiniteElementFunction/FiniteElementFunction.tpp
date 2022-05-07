/*
 * FiniteElementFunction.tpp
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

#ifndef FINITE_ELEMENT_FUNCTION_TPP
#define FINITE_ELEMENT_FUNCTION_TPP

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.hpp>
#include <LightFEM/Expression/Function/Const.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunction.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

template<ExprType Type> template<typename Expr>
FiniteElementFunction<Type>::FiniteElementFunction(const FiniteElementFunctionExpression< Type, Expr >& expr) :
	m_fSpace(expr.getFunctionSpace())
{
	m_isDiscretized = false;
	m_coefs.resize(m_fSpace->getNBasisFunction());

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }
}

template<ExprType Type> template<typename Expr>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::operator=  (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	m_fSpace = expr.getFunctionSpace();
	m_coefs.resize(m_fSpace->getNBasisFunction());

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::operator+= (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] += expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::operator-= (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] -= expr.getCoef(i); }

	return *this;
}

template<ExprType Type>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::operator*= (const Const& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] *= expr.getValue(); }

	return *this;
}

template<ExprType Type>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::operator/= (const Const& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] /= expr.getValue(); }

	return *this;
}

template<ExprType Type>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::discretize()
{
	m_values.resize(getMesh()->getNElem());
	m_grad.resize(getMesh()->getNElem());

	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		m_values[e].setElement(getMesh()->getElem(e));
		m_grad[e].setElement(getMesh()->getElem(e));

		const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(-1);
		const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(-1);
		for (size_t k=0;k<Element::getNxiNd();++k)
		{
			m_values[e][k] = values[k];
			m_grad[e][k] = grad[k];
		}
	}

#pragma omp parallel for
	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		for (size_t locId=0;locId<m_fSpace->getNBasisFunctionPerElement();++locId)
		{
			const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(locId);
			const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(locId);

			for (size_t k=0;k<Element::getNxiNd();++k)
			{
				m_values[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*values[k];
				m_grad[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*(transpose(getMesh()->getElem(e)->getInvJacobianDisc(k))*grad[k]);
			}
		}
	}
	m_isDiscretized = true;

	return *this;
}

template<ExprType Type>
FiniteElementFunction<Type>& FiniteElementFunction<Type>::discretize(MPI_Comm com)
{
	m_values.resize(getMesh()->getNElem());
	m_grad.resize(getMesh()->getNElem());

	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		m_values[e].setElement(getMesh()->getElem(e));
		m_grad[e].setElement(getMesh()->getElem(e));

		const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(-1);
		const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(-1);
		for (size_t k=0;k<Element::getNxiNd();++k)
		{
			m_values[e][k] = values[k];
			m_grad[e][k] = grad[k];
		}
	}

	#pragma omp parallel
	{
		MpiRange range(com, getMesh()->getNElem());
		
		#pragma omp for
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t locId=0;locId<m_fSpace->getNBasisFunctionPerElement();++locId)
			{
				const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(locId);
				const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(locId);

				for (size_t k=0;k<Element::getNxiNd();++k)
				{
					m_values[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*values[k];
					m_grad[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*(transpose(getMesh()->getElem(e)->getInvJacobianDisc(k))*grad[k]);
				}
			}
		}
	}
	m_isDiscretized = true;

	return *this;
}

////////////////////////////////////////////////////////////////////////

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>::CpxFiniteElementFunction(const FiniteElementFunctionExpression< Type, Expr >& expr) :
	m_fSpace(expr.getFunctionSpace())
{
	m_isDiscretized = false;
	m_coefs.resize(m_fSpace->getNBasisFunction());

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }
}


template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>::CpxFiniteElementFunction(const CpxFiniteElementFunctionExpression< Type, Expr >& expr) :
	m_fSpace(expr.getFunctionSpace())
{
	m_isDiscretized = false;
	m_coefs.resize(m_fSpace->getNBasisFunction());

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator=  (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	m_fSpace = expr.getFunctionSpace();

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator=  (const CpxFiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	m_fSpace = expr.getFunctionSpace();

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] = expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator+= (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] += expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator+= (const CpxFiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] += expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator-= (const FiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] -= expr.getCoef(i); }

	return *this;
}

template<ExprType Type> template<typename Expr>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator-= (const CpxFiniteElementFunctionExpression< Type, Expr >& expr)
{
	m_isDiscretized = false;
	if (m_fSpace != expr.getFunctionSpace()) { throw std::invalid_argument("Both functions must be in the the same function space."); }

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] -= expr.getCoef(i); }

	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator*= (const Const& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] *= expr.getValue(); }

	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator*= (const CpxConst& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] *= expr.getValue(); }

	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator/= (const Const& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] /= expr.getValue(); }

	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::operator/= (const CpxConst& expr)
{
	m_isDiscretized = false;

	for (size_t i=0;i<m_fSpace->getNBasisFunction();++i) { m_coefs[i] /= expr.getValue(); }

	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::discretize()
{
	m_values.resize(getMesh()->getNElem());
	m_grad.resize(getMesh()->getNElem());

	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		m_values[e].setElement(getMesh()->getElem(e));
		m_grad[e].setElement(getMesh()->getElem(e));

		const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(-1);
		const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(-1);
		for (size_t k=0;k<Element::getNxiNd();++k)
		{
			m_values[e][k] = values[k];
			m_grad[e][k] = grad[k];
		}
	}

#pragma omp parallel for
	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		for (size_t locId=0;locId<m_fSpace->getNBasisFunctionPerElement();++locId)
		{
			const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(locId);
			const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(locId);

			for (size_t k=0;k<Element::getNxiNd();++k)
			{
				m_values[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*values[k];
				m_grad[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*(transpose(getMesh()->getElem(e)->getInvJacobian(k))*grad[k]);
			}
		}
	}
	m_isDiscretized = true;
	return *this;
}

template<ExprType Type>
CpxFiniteElementFunction<Type>& CpxFiniteElementFunction<Type>::discretize(MPI_Comm com)
{
	m_values.resize(getMesh()->getNElem());
	m_grad.resize(getMesh()->getNElem());

	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		m_values[e].setElement(getMesh()->getElem(e));
		m_grad[e].setElement(getMesh()->getElem(e));

		const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(-1);
		const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(-1);
		for (size_t k=0;k<Element::getNxiNd();++k)
		{
			m_values[e][k] = values[k];
			m_grad[e][k] = grad[k];
		}
	}

	#pragma omp parallel
	{
		MpiRange range(com, getMesh()->getNElem());
		
		#pragma omp for
		for (size_t e=range.begin();e<range.end();++e)
		{
			for (size_t locId=0;locId<m_fSpace->getNBasisFunctionPerElement();++locId)
			{
				const std::vector< typename FSpace<Type>::ValueType >& values = m_fSpace->getDiscBasisFunction(locId);
				const std::vector< typename FSpace<Type>::GradValueType >& grad = m_fSpace->getDiscGradXiBasisFunction(locId);

				for (size_t k=0;k<Element::getNxiNd();++k)
				{
					m_values[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*values[k];
					m_grad[e][k] += m_coefs[m_fSpace->getGlobalId(e, locId)]*(transpose(getMesh()->getElem(e)->getInvJacobian(k))*grad[k]);
				}
			}
		}
	}
	m_isDiscretized = true;
	return *this;
}

#endif // FINITE_ELEMENT_FUNCTION_TPP
