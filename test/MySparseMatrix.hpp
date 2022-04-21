#ifndef MY_SPARSE_MATRIX_HPP
#define MY_SPARSE_MATRIX_HPP

#include <vector>
#include <utility>
#include <tuple>

template<typename T>
class MySparseMatrix
{
public:
	typedef std::tuple<int, int, T> Entry;
	typedef std::vector< std::pair< int, T > > RowSlice;
public:
	MySparseMatrix(const int nrows, const int ncols) : m_rowSlice(nrows), m_ncols(ncols) {}
	void setFromTriplets(typename std::vector< Entry >::const_iterator begin_entries, typename std::vector< Entry >::const_iterator end_entries);
	
	std::vector< T > apply(const std::vector< T >& x) const; // returns Mx
	T  apply(const std::vector< T >& x, int i) const; // returns Mx_i
	
	int getNrows() const { return m_rowSlice.size(); }
	int getNCols() const { return m_ncols; }
private:
	std::vector< RowSlice > m_rowSlice;
	int m_ncols;
};

template<typename T>
void MySparseMatrix<T>::setFromTriplets(typename std::vector< Entry >::const_iterator begin_entries, typename std::vector< Entry >::const_iterator end_entries)
{
	for (typename std::vector< Entry >::const_iterator entryIt = begin_entries; entryIt != end_entries; ++entryIt)
	{
		const auto& [i, j, aij] = *entryIt;
		//std::cout << "A(" << i << ", " << j << ") = " << aij << std::endl;
		
		//std::cout << "m_rowSlice[" << i << "].size() = " << m_rowSlice[i].size() << std::endl;
		
		RowSlice& Ai = m_rowSlice[i];
		//std::cout << "A[" << i << ",:] = {" << std::flush;
		//for (size_t idx=0;idx<Ai.size();++idx) { std::cout << "(" << Ai[idx].first << ", " << Ai[idx].second << ")" << std::endl; } 
		//std::cout << "}" << std::endl;
		typename RowSlice::const_iterator it = std::find_if(std::cbegin(Ai), std::cend(Ai), [j](const std::pair< int, T >& pair)
		{
			return pair.first == j;
		});
		if (it == std::cend(Ai)) { Ai.push_back(std::make_pair(j, aij)); }
		else { Ai[it->first].second += aij; }
	}
	for (size_t i=0;i<m_rowSlice.size();++i)
	{
		m_rowSlice[i].shrink_to_fit();
	} 
}

template<typename T>
std::vector< T > MySparseMatrix<T>::apply(const std::vector< T >& x) const
{
	std::vector< T > Mx(getNrows());
	
	if (x.size() != getNCols()) { throw std::invalid_argument("Matrix and vector dimension must match."); }
	
	#pragma omp parallel for
	for (int i=0;i<m_rowSlice.size();++i)
	{
		for (const auto& [j, aij] : m_rowSlice[i])
		{
			Mx[i] += aij*x[j];
		}
	}
	
	return Mx;
}

template<typename T>
T MySparseMatrix<T>::apply(const std::vector< T >& x, int i) const
{
	T Mxi = 0.0;
	
	if (int(x.size()) != getNCols()) { throw std::invalid_argument("Matrix and vector dimension must match."); }
	

	for (const auto& [j, aij] : m_rowSlice[i])
	{
		Mxi += aij*x[j];
	}
	
	return Mxi;
}

#endif // MY_SPARSE_MATRIX_HPP
