#include <LightFEM/Expression/LinAlg/Matrix.hpp>

Matrix::Matrix(std::initializer_list< std::initializer_list< double > > values) : 
	m_nrows(values.size()),
	m_ncols(std::cbegin(values)->size()),
	m_core(m_nrows*m_ncols)
{
	size_t i=0;
	for (const std::initializer_list< double >& row_values : values)
	{
		size_t j=0;
		if (row_values.size() != m_ncols) { throw std::invalid_argument("Row " + std::to_string(i) + " doesn't have as many values as the first row."); }
		for (const double& value: row_values)
		{
			m_core[flatIndex(i,j)] = value;
			++j;
		}
		++i;
	}
}

////////////////////////////////////////////////////////////////////////

CpxMatrix::CpxMatrix(std::initializer_list< std::initializer_list< std::complex< double > > > values) : 
	m_nrows(values.size()),
	m_ncols(std::cbegin(values)->size()),
	m_core(m_nrows*m_ncols)
{
	size_t i=0;
	for (const std::initializer_list< std::complex< double > >& row_values : values)
	{
		size_t j=0;
		if (row_values.size() != m_ncols) { throw std::invalid_argument("Row " + std::to_string(i) + " doesn't have as many values as the first row."); }
		for (const std::complex< double >& value: row_values)
		{
			m_core[flatIndex(i,j)] = value;
			++j;
		}
		++i;
	}
}
