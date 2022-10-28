#include <LightFEM/Mesh/operators.hpp>

void printFunctionBin(const std::string& fname, const ScalarField& expr)
{
	std::ofstream out(fname, std::ios_base::binary);
	
	for (std::size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		out.write( reinterpret_cast<const char *>(expr[e].data()), Element::getNxiNd()*sizeof(ElementWiseScalarField::ValueType) );
	}
}

void printFunctionBin(const std::string& fname, const CpxScalarField& expr)
{
	std::ofstream out(fname, std::ios_base::binary);
	
	for (std::size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		out.write( reinterpret_cast<const char *>(expr[e].data()), Element::getNxiNd()*sizeof(CpxElementWiseScalarField::ValueType) );
	}
}

ScalarField readFunctionBin(const std::string& fname, const Mesh* mesh)
{
	ScalarField expr(mesh);
	
	std::ifstream in(fname, std::ios_base::binary);

	if (!in.is_open()) { throw std::runtime_error("failed to open " + fname); }
	
	for (std::size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		in.read( reinterpret_cast<char *>(expr[e].data()), Element::getNxiNd()*sizeof(ElementWiseScalarField::ValueType) );
	}
	
	return expr;
}

CpxScalarField readCpxFunctionBin(const std::string& fname, const Mesh* mesh)
{
	CpxScalarField expr(mesh);
	
	std::ifstream in(fname, std::ios_base::binary);

	if (!in.is_open()) { throw std::runtime_error("failed to open " + fname); }
	
	for (std::size_t e=0;e<expr.getMesh()->getNElem();++e)
	{
		in.read( reinterpret_cast<char *>(expr[e].data()), Element::getNxiNd()*sizeof(CpxElementWiseScalarField::ValueType) );
	}
	
	return expr;
}
