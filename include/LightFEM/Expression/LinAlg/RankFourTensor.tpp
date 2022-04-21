/*
 * RankFourTensor.tpp
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

#ifndef RANK_FOUR_TENSOR_TPP
#define RANK_FOUR_TENSOR_TPP

#include <LightFEM/Expression/LinAlg/RankFourTensor.hpp>

#include <numeric>
#include <functional>

template<typename Expr>
RankFourTensor::RankFourTensor(const RankFourTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
}

template<typename Expr>
RankFourTensor& RankFourTensor::operator= (const RankFourTensorExpression< Expr >& expr)
{ 
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankFourTensor& RankFourTensor::operator+= (const RankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] +=  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankFourTensor& RankFourTensor::operator-= (const RankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] -= expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankFourTensor& RankFourTensor::operator*= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] *=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankFourTensor& RankFourTensor::operator/= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] /=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

////////////////////////////////////////////////////////////////////////

template<typename Expr>
CpxRankFourTensor::CpxRankFourTensor(const RankFourTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
}

template<typename Expr>
CpxRankFourTensor::CpxRankFourTensor(const CpxRankFourTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator= (const RankFourTensorExpression< Expr >& expr)
{
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator= (const CpxRankFourTensorExpression< Expr >& expr)
{
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] =  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator+= (const RankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] +=  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator+= (const CpxRankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] +=  expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator-= (const RankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] -= expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator-= (const CpxRankFourTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] -= expr(i,j,k,l); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator*= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] *=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator*= (const CpxScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] *=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator/= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] /=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankFourTensor& CpxRankFourTensor::operator/= (const CpxScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }
	
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			for (size_t k=0;k<m_shape[2];++k) 
			{ 
				for (size_t l=0;l<m_shape[3];++l) 
				{ 
					m_core[flatIndex(i,j,k,l)] /=  expr.eval(); 
				}
			}
		} 
	} 
	
	return *this;
}

#endif // RANK_FOUR_TENSOR_TPP
