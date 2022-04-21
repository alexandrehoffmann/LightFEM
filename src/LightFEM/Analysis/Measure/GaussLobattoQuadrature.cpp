/*
 * GaussLobattoQuadrature.cpp
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

#include <LightFEM/Analysis/Measure/GaussLobattoQuadrature.hpp>

GaussLobattoQuadrature::GaussLobattoQuadrature(const Mesh *mesh, const size_t N) :
	Measure(mesh)
{
	std::valarray< long double > xi = getNodes(N);
	std::valarray< long double > wi(N);

	for (size_t i=0;i<N;++i)
	{
		const long double Pnm1 = P(N-1, xi[i]);
		wi[i] = 2.0 / (N*(N-1)*(Pnm1*Pnm1));
	}

	for (std::size_t e=0;e<m_mesh->getNElem();++e)
	{
		m_nodesAndWeightsPerElement[e].resize(N*N);
		m_boundaryNodesAndWeightsPerElement[e].resize(N);
		for (size_t i=0;i<N;++i) { for (size_t j=0;j<N;++j)
		{
			m_nodesAndWeightsPerElement[e][i*N + j].xi1 = xi[i];
			m_nodesAndWeightsPerElement[e][i*N + j].xi2 = xi[j];
			m_nodesAndWeightsPerElement[e][i*N + j].wi  = wi[i]*wi[j]*m_mesh->getElem(e)->getAbsDetJacobian(xi[i], xi[j]);
		}}
		for (size_t i=0;i<N;++i)
		{
			m_boundaryNodesAndWeightsPerElement[e][i].t        = xi[i];
			m_boundaryNodesAndWeightsPerElement[e][i].topWi    = wi[i]*m_mesh->getElem(e)->getDs(Element::Boundary::TOP    , xi[i]);
			m_boundaryNodesAndWeightsPerElement[e][i].bottomWi = wi[i]*m_mesh->getElem(e)->getDs(Element::Boundary::BOTTTOM, xi[i]);
			m_boundaryNodesAndWeightsPerElement[e][i].leftWi   = wi[i]*m_mesh->getElem(e)->getDs(Element::Boundary::LEFT   , xi[i]);
			m_boundaryNodesAndWeightsPerElement[e][i].rightWi  = wi[i]*m_mesh->getElem(e)->getDs(Element::Boundary::RIGHT  , xi[i]);
		}
	}
}

std::valarray< long double > GaussLobattoQuadrature::getNodes(const size_t N) //c.f. Aberth method
{
	constexpr long double epsilon = 1.0e-14;
	auto eq = [epsilon](const long double a, const long double b) -> bool { return (fabs(a - b) <= epsilon * std::max(1.0L, std::max(a, b))); };
	auto leq = [&eq](const long double a, const long double b) -> bool { return (a < b) or eq(a,b); };
	auto geq = [&eq](const long double a, const long double b) -> bool { return (a > b) or eq(a,b); };

	std::valarray< long double > err(N);
	std::valarray< long double > xi(N);
	for (std::size_t i=0;i<N;++i)
	{
		xi[i] = -1.0 + 2.0*double(i) / double(N-1);
		err[i] = fabs(dP(N, xi[i]));
	}

	while (not eq(err.max(), 0.0))
	{
		std::valarray< long double > dxi(N);
		std::valarray< long double > alpha(1.0, N);
		for (size_t i=1;i<N-1;++i)
		{
			const long double rho = dP(N-1, xi[i]) / d2P(N-1, xi[i]);
			long double sum = 0.0;
			for (std::size_t j=0;j<N;++j) { if (j != i)
			{
				sum += 1.0 / (xi[i] - xi[j]);
			}}
			dxi[i] = -rho / (1.0 - rho*sum);
			if (leq(xi[i] + dxi[i], -1.0)) { alpha[i] = (-1.0 + 1.0e-3 - xi[i]) / dxi[i]; }
			if (geq(xi[i] + dxi[i],  1.0)) { alpha[i] = ( 1.0 - 1.0e-3 - xi[i]) / dxi[i]; }
		}
		const long double min_alpha = alpha.min();
		for (size_t i=0;i<N;++i)
		{
			xi[i] += min_alpha*dxi[i];
			err[i] = fabs(dP(N-1, xi[i]));
		}
	}

	std::sort(std::begin(xi), std::end(xi));

	return xi;
}

long double GaussLobattoQuadrature::dP(const size_t N, const long double x)
{
	constexpr long double epsilon = 1.0e-14;
	auto eq = [epsilon](const long double a, const long double b) -> bool { return (fabs(a - b) <= epsilon * std::max(1.0L, std::max(a, b))); };

	if (eq(fabs(x), 1.0)) { return 0.0; }
	else if (N == 0) { return 1.0; }

	return ( N*P(N-1,x) - N*x*P(N,x) ) / (1.0 - x*x); // c.f. https://mathworld.wolfram.com/LegendrePolynomial.html
}

long double GaussLobattoQuadrature::d2P(const size_t N, const long double x)
{
	constexpr long double epsilon = 1.0e-14;
	auto eq = [epsilon](const long double a, const long double b) -> bool { return (fabs(a - b) <= epsilon * std::max(1.0L, std::max(a, b))); };

	if (N == 0) { return 0.0; }
	if (eq(fabs(x), 1.0)) { return 1.0; }
//	return (2.0*x*dP(N, x) - long double(N*(N-1))*P(N,x)) / ((1.0 - x*x)*(1.0 - x*x));
	const long double u = N*P(N-1,x) - N*x*P(N,x);
	const long double du = N*dP(N-1,x) - N*x*dP(N,x) - N*P(N,x);
	const long double v = 1.0 - x*x;
	const long double dv = -2.0*x;

	return (du*v - dv*u) / (v*v);
}
