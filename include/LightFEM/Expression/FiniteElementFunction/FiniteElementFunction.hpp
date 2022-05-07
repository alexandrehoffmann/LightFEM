/*
 * FiniteElementFunction.hpp
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

#ifndef FINITE_ELEMENT_FUNCTION_HPP
#define FINITE_ELEMENT_FUNCTION_HPP

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Tools/MpiRange.hpp>

class functionNotDiscretized : public std::exception
{
	virtual const char* what() const throw()
	{
		return "function must be discretized first";
	}
};

template<ExprType Type>
class FiniteElementFunction : public FiniteElementFunctionExpression< Type, FiniteElementFunction<Type> >
{
public:
	typedef typename Traits< FiniteElementFunction<Type> >::ReturnType ReturnType;
	typedef typename Traits< FiniteElementFunction<Type> >::GradReturnType GradReturnType;
	typedef typename Traits< FiniteElementFunction<Type> >::DiffReturnType DiffReturnType;

	typedef typename Traits< FiniteElementFunction<Type> >::ValueType ValueType;
	typedef typename Traits< FiniteElementFunction<Type> >::GradValueType GradValueType;

	typedef typename FSpace< Type >::Type FSpaceType;
public:
	FiniteElementFunction(const FSpaceType* fSpace) : m_fSpace(fSpace), m_coefs(fSpace->getNBasisFunction()) {}
	FiniteElementFunction(const FSpaceType* fSpace, const Scalar& coef) : m_fSpace(fSpace), m_coefs(fSpace->getNBasisFunction(), coef) {}
	FiniteElementFunction(const FSpaceType* fSpace, const Scalar* coefs) : m_fSpace(fSpace), m_coefs(coefs, coefs+fSpace->getNBasisFunction()) {}
	template<typename Expr> FiniteElementFunction(const FiniteElementFunctionExpression< Type, Expr >& expr);

	template<typename Expr> FiniteElementFunction<Type>& operator=  (const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> FiniteElementFunction<Type>& operator+= (const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> FiniteElementFunction<Type>& operator-= (const FiniteElementFunctionExpression< Type, Expr >& expr);
	FiniteElementFunction<Type>& operator*= (const Const& expr);
	FiniteElementFunction<Type>& operator/= (const Const& expr);
public:
	inline ReturnType     operator[] (const size_t e)                   const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return m_values[e]; }
	inline GradReturnType getGrad    (const size_t e)                   const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return m_grad[e]; }
	inline DiffReturnType getD       (const size_t dim, const size_t e) const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return DiffReturnType(getGrad(e), dim); }
	inline const Scalar&  getCoef    (const size_t globId)              const { return m_coefs[globId]; }
	
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }

	inline void setCoef(const size_t globId, const Scalar& c) { m_isDiscretized = false; m_coefs[globId] = c; }

	inline const FSpaceType* getFunctionSpace() const { return m_fSpace; }
	inline const Mesh*       getMesh()          const { return m_fSpace->getMesh(); }

	FiniteElementFunction<Type>& discretize();
	/**
	 * The discretizeation isn't broadcasted over MPI processes. If you want to use this in a bilinear or linear form please use the same MPI_Comm
	 * if you want to compute the integral of this function or an expression involving this function, please use MPI_Allreduce with the same MPI_Comm
	 */
	FiniteElementFunction<Type>& discretize(MPI_Comm com); 
	
	inline bool isDiscretized() const { return m_isDiscretized; } 
private:
	std::vector< ValueType >     m_values;
	std::vector< GradValueType > m_grad;
	bool m_isDiscretized;
	
	const FSpaceType* m_fSpace;
	std::vector< Scalar > m_coefs;
};

typedef FiniteElementFunction<ExprType::SCALAR> FiniteElementScalarField;
typedef FiniteElementFunction<ExprType::VECTOR> FiniteElementVectorField;

////////////////////////////////////////////////////////////////////////

template<ExprType Type>
class CpxFiniteElementFunction : public CpxFiniteElementFunctionExpression< Type, CpxFiniteElementFunction<Type> >
{
public:
	typedef typename Traits< CpxFiniteElementFunction<Type> >::ReturnType ReturnType;
	typedef typename Traits< CpxFiniteElementFunction<Type> >::GradReturnType GradReturnType;
	typedef typename Traits< CpxFiniteElementFunction<Type> >::DiffReturnType DiffReturnType;

	typedef typename Traits< CpxFiniteElementFunction<Type> >::ValueType ValueType;
	typedef typename Traits< CpxFiniteElementFunction<Type> >::GradValueType GradValueType;

	typedef typename FSpace< Type >::Type FSpaceType;
public:
	CpxFiniteElementFunction(const FSpaceType* fSpace) : m_fSpace(fSpace), m_coefs(fSpace->getNBasisFunction()) {}
	CpxFiniteElementFunction(const FSpaceType* fSpace, const Scalar& coef) : m_fSpace(fSpace), m_coefs(fSpace->getNBasisFunction(), coef) {}
	CpxFiniteElementFunction(const FSpaceType* fSpace, const CpxScalar* coefs) : m_fSpace(fSpace), m_coefs(coefs, coefs+fSpace->getNBasisFunction()) {}
	template<typename Expr> CpxFiniteElementFunction(const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction(const CpxFiniteElementFunctionExpression< Type, Expr >& expr);

	template<typename Expr> CpxFiniteElementFunction<Type>& operator=  (const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction<Type>& operator=  (const CpxFiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction<Type>& operator+= (const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction<Type>& operator+= (const CpxFiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction<Type>& operator-= (const FiniteElementFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFiniteElementFunction<Type>& operator-= (const CpxFiniteElementFunctionExpression< Type, Expr >& expr);
	CpxFiniteElementFunction<Type>& operator*= (const Const& expr);
	CpxFiniteElementFunction<Type>& operator*= (const CpxConst& expr);
	CpxFiniteElementFunction<Type>& operator/= (const Const& expr);
	CpxFiniteElementFunction<Type>& operator/= (const CpxConst& expr);
public:
	inline ReturnType       operator[] (const size_t e)                   const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return m_values[e]; }
	inline GradReturnType   getGrad    (const size_t e)                   const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return m_grad[e]; }
	inline DiffReturnType   getD       (const size_t dim, const size_t e) const { if (not m_isDiscretized) { throw functionNotDiscretized(); } return DiffReturnType(getGrad(e), dim); }
	inline const CpxScalar& getCoef    (const size_t globId)              const { return m_coefs[globId]; }
	
	inline void setCoef(const size_t globId, const Scalar& c) { m_isDiscretized = false; m_coefs[globId] = c; }
	inline void setCoef(const size_t globId, const CpxScalar& c) { m_isDiscretized = false; m_coefs[globId] = c; }

	inline const FSpaceType* getFunctionSpace() const { return m_fSpace; }
	inline const Mesh*       getMesh()          const { return m_fSpace->getMesh(); }

	CpxFiniteElementFunction<Type>& discretize();
	/**
	 * The discretizeation isn't broadcasted over MPI processes. If you want to use this in a bilinear or linear form please use the same MPI_Comm
	 * if you want to compute the integral of this function or an expression involving this function, please use MPI_Allreduce with the same MPI_Comm
	 */
	CpxFiniteElementFunction<Type>& discretize(MPI_Comm com); 
	inline bool isDiscretized() const { return m_isDiscretized; } 
private:
	std::vector< ValueType >     m_values;
	std::vector< GradValueType > m_grad;
	bool m_isDiscretized;
	
	const FSpaceType* m_fSpace;
	std::vector< CpxScalar > m_coefs;
};

typedef CpxFiniteElementFunction<ExprType::SCALAR> CpxFiniteElementScalarField;
typedef CpxFiniteElementFunction<ExprType::VECTOR> CpxFiniteElementVectorField;

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.tpp>

#endif // FINITE_ELEMENT_FUNCTION_HPP
