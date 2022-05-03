/*
 * Element.cpp
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

#include <LightFEM/Mesh/Element.hpp>

#include <LightFEM/Expression/LinAlg/operators.hpp>

#include <LightFEM/Tools/FixedSizeMatrix.hpp>
#include <LightFEM/Tools/FixedSizeVector.hpp>
#include <LightFEM/Tools/FixedSizeLuSolver.hpp>

const std::array<std::array<NodeRef, Element::ORDER+1>, Element::ORDER+1> Element::Xi = {
	std::array<NodeRef, Element::ORDER+1>{NodeRef(-1.0, -1.0), NodeRef(-1.0, 0.0), NodeRef(-1.0, 1.0)},
	std::array<NodeRef, Element::ORDER+1>{NodeRef( 0.0, -1.0), NodeRef( 0.0, 0.0), NodeRef( 0.0, 1.0)},
	std::array<NodeRef, Element::ORDER+1>{NodeRef( 1.0, -1.0), NodeRef( 1.0, 0.0), NodeRef( 1.0, 1.0)}
};

Element::Element(const std::array<std::array<NodeWorld*, 2>, 2>& X) :
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4)
{
	init(X);
}

Element::Element(const std::array<std::array<NodeWorld*, 2>, 2>& X, std::initializer_list< int > domainIds) :
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4),
	m_domainIds(domainIds)
{
	init(X);
}

Element::Element(const std::array<std::array<NodeWorld*, 2>, 2>& X, const std::vector< int >& domainIds) :
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4),
	m_domainIds(domainIds)
{
	init(X);
}

void Element::init(const std::array<std::array<NodeWorld*, 2>, 2>& X)
{
	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		m_deleteX[i][j] = false;
	}}

	m_X[0][1] = new NodeWorld(); m_deleteX[0][1] = true;
	m_X[1][0] = new NodeWorld(); m_deleteX[1][0] = true;
	m_X[1][1] = new NodeWorld(); m_deleteX[1][1] = true;
	m_X[1][2] = new NodeWorld(); m_deleteX[1][2] = true;
	m_X[2][1] = new NodeWorld(); m_deleteX[2][1] = true;

	m_X[0][0] = X[0][0];
	m_X[0][2] = X[0][1];
	m_X[0][1]->x = 0.5*(m_X[0][0]->x + m_X[0][2]->x);
	m_X[0][1]->y = 0.5*(m_X[0][0]->y + m_X[0][2]->y);

	m_X[2][0] = X[1][0];
	m_X[2][2] = X[1][1];
	m_X[2][1]->x = 0.5*(m_X[2][0]->x + m_X[2][2]->x);
	m_X[2][1]->y = 0.5*(m_X[2][0]->y + m_X[2][2]->y);

	m_X[1][0]->x = 0.5*(m_X[0][0]->x + m_X[2][0]->x);
	m_X[1][0]->y = 0.5*(m_X[0][0]->y + m_X[2][0]->y);
	m_X[1][1]->x = 0.5*(m_X[0][1]->x + m_X[2][1]->x);
	m_X[1][1]->y = 0.5*(m_X[0][1]->y + m_X[2][1]->y);
	m_X[1][2]->x = 0.5*(m_X[0][2]->x + m_X[2][2]->x);
	m_X[1][2]->y = 0.5*(m_X[0][2]->y + m_X[2][2]->y);

	initTransport();
	for (size_t i=0;i<xi.size();++i) { for (size_t j=0;j<xi.size();++j)
	{
		m_abs_dejJ[index2d(i,j)] = getAbsDetJacobian(xi[i], xi[j]);
		m_invJ[index2d(i,j)] = getInvJacobian(xi[i], xi[j]);
	}}

	for (const Boundary& boundary : AllBoundaries) { for (size_t k=0;k<xi.size();++k)
	{
		m_normal[index2d(boundary, k)] = getNormal(boundary, xi[k]);
		m_ds[index2d(boundary, k)] = getDs(boundary, xi[k]);
	}}
}

Element::Element(const std::array<std::array<NodeWorld*, 3>, 3>& X) :
	m_X(X),
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4)
{
	init();
}

Element::Element(const std::array<std::array<NodeWorld*, 3>, 3>& X, std::initializer_list< int > domainIds) :
	m_X(X),
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4),
	m_domainIds(domainIds)
{
	init();
}

Element::Element(const std::array<std::array<NodeWorld*, 3>, 3>& X, const std::vector< int >& domainIds) :
	m_X(X),
	m_abs_dejJ(getNxiNd()),
	m_ds(getNxi()*4),
	m_invJ(getNxiNd()),
	m_normal(getNxi()*4),
	m_domainIds(domainIds)
{
	init();
}

void Element::init()
{
	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		m_deleteX[i][j] = false;
	}}

	initTransport();
	for (size_t i=0;i<xi.size();++i) { for (size_t j=0;j<xi.size();++j)
	{
		m_abs_dejJ[index2d(i,j)] = getAbsDetJacobian(xi[i], xi[j]);
		m_invJ[index2d(i,j)] = getInvJacobian(xi[i], xi[j]);
	}}

	for (const Boundary& boundary : AllBoundaries) { for (size_t k=0;k<xi.size();++k)
	{
		m_normal[index2d(boundary, k)] = getNormal(boundary, xi[k]);
		m_ds[index2d(boundary, k)] = getDs(boundary, xi[k]);
	}}
}

Element::~Element()
{
	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
			if (m_deleteX[i][j]) { delete m_X[i][j]; }
	}}
}

NodeWorld Element::getXworld(const NodeRef& Xi) const
{
	NodeWorld X(0.0, 0.0);
	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		X.x += m_c[i][j].x*std::pow(Xi.xi1, i)*std::pow(Xi.xi2, j);
		X.y += m_c[i][j].y*std::pow(Xi.xi1, i)*std::pow(Xi.xi2, j);
	}}
	return X;
}

Matrix Element::getJacobian(const NodeRef& Xi) const
{
	Matrix J(2,2);
	// derivatives with respect to xi1
	auto [xi1, xi2] = Xi;

	for (size_t i=1;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		J(0,0) += double(i)*m_c[i][j].x*std::pow(xi1, i-1)*std::pow(xi2, j);
		J(1,0) += double(i)*m_c[i][j].y*std::pow(xi1, i-1)*std::pow(xi2, j);
	}}
	// derivatives with respect to xi2
	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=1;j<ORDER+1;++j)
	{
		J(0,1) += double(j)*m_c[i][j].x*std::pow(xi1, i)*std::pow(xi2, j-1);
		J(1,1) += double(j)*m_c[i][j].y*std::pow(xi1, i)*std::pow(xi2, j-1);
	}}

	return J;
}

Vector Element::getNormal(const Boundary b, const double t) const
{
	std::function<NodeRef(double)> t_to_xi;
	size_t colId;
	
	Matrix R(2,2);
	
	if (b == Boundary::TOP)
	{
		t_to_xi = [](const double t) -> NodeRef { return NodeRef(t, 1.0); };
		colId = 0;
		R(0,0) = 0.0; R(0,1) = -1.0;
		R(1,0) = 1.0; R(1,1) = 0.0;
	}
	else if (b == Boundary::BOTTTOM)
	{
		t_to_xi = [](const double t) -> NodeRef { return NodeRef(t, -1.0); };
		colId = 0;
		R(0,0) = 0.0; R(0,1) = 1.0;
		R(1,0) = -1.0; R(1,1) = 0.0;
	}
	else if (b == Boundary::LEFT)
	{
		t_to_xi = [](const double t) -> NodeRef { return NodeRef(-1.0, t); };
		colId = 1;
		R(0,0) = 0.0; R(0,1) = -1.0;
		R(1,0) = 1.0; R(1,1) = 0.0;
	}
	else
	{
		t_to_xi = [](const double t) -> NodeRef { return NodeRef(1.0, t); };
		colId = 1;
		R(0,0) = 0.0; R(0,1) = 1.0;
		R(1,0) = -1.0; R(1,1) = 0.0;
	}
	
	const Matrix J = getJacobian(t_to_xi(t));
	
	const Vector tan{J(0,colId), J(1,colId)};
	return R*tan / norm(tan);
}

double Element::getAbsDetJacobian(const NodeRef& Xi) const
{
	Matrix J = getJacobian(Xi);
	return std::fabs(J(0,0)*J(1,1) - J(0,1)* J(1,0));
}

double Element::getDs(const Boundary b, const double t) const
{
	std::function<NodeRef(double)> t_to_xi;
	size_t colId;
	
	if (b == Boundary::TOP)          { t_to_xi = [](const double t) -> NodeRef { return NodeRef(t, 1.0); } ; colId = 0; }
	else if (b == Boundary::BOTTTOM) { t_to_xi = [](const double t) -> NodeRef { return NodeRef(t, -1.0); }; colId = 0; }
	else if (b == Boundary::LEFT)    { t_to_xi = [](const double t) -> NodeRef { return NodeRef(-1.0, t); }; colId = 1; }
	else if (b == Boundary::RIGHT) 	 { t_to_xi = [](const double t) -> NodeRef { return NodeRef(1.0, t); };  colId = 1; }
	
	const Matrix J = getJacobian(t_to_xi(t));
	const Vector tan{J(0,colId), J(1,colId)};

	return norm(tan);
}

Matrix Element::getInvJacobian(const NodeRef& Xi) const
{
	Matrix J = getJacobian(Xi);

	const double detJ = J(0,0)*J(1,1) - J(0,1)* J(1,0);
	Matrix invJ(2,2);

	invJ(0,0) =  J(1,1) / detJ;
	invJ(0,1) = -J(0,1) / detJ;
	invJ(1,0) = -J(1,0) / detJ;
	invJ(1,1) =  J(0,0) / detJ;

	return invJ;
}

void Element::initTransport()
{
	// T_x(Xi) = sum_kl c_kl1 Xi_1^k Xi_2^l
	// T_y(Xi) = sum_kl c_kl2 Xi_1^k Xi_2^l
	// we must solve T_x(Xi_ij) = sum_kl Xi_ij_1^k Xi_ij_2^l c_kl_x = X_ij_x
	// we must solve T_y(Xi_ij) = sum_kl Xi_ij_1^k Xi_ij_2^l c_kl_y = X_ij_y
	constexpr size_t N = (ORDER+1)*(ORDER+1);

	FixedSizeMatrix<N,N> M;
	FixedSizeVector<N> x_world;
	FixedSizeVector<N> y_world;

	auto index2d = [](size_t i, size_t j) -> size_t { return i*(ORDER+1) + j; };

	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		x_world[index2d(i,j)] = m_X[i][j]->x;
		y_world[index2d(i,j)] = m_X[i][j]->y;

		for (size_t k=0;k<ORDER+1;++k) { for (size_t l=0;l<ORDER+1;++l)
		{
			M(index2d(i,j), index2d(k,l)) = std::pow(Xi[i][j].xi1, k)*std::pow(Xi[i][j].xi2, l);
		}}
	}}
	FixedSizeLuSolver<N> LU(M);
	FixedSizeVector<N> c_x = LU.solve(x_world);
	FixedSizeVector<N> c_y = LU.solve(y_world);

	for (size_t i=0;i<ORDER+1;++i) { for (size_t j=0;j<ORDER+1;++j)
	{
		m_c[i][j].x = c_x[index2d(i,j)];
		m_c[i][j].y = c_y[index2d(i,j)];
	}}
}

#include <LightFEM/Tools/Linesearch.hpp>

std::pair< bool, NodeRef > Element::getXRef(const NodeWorld& Xworld, const double epsilon) const
{
	// we want to find Xi such that 
	// T_x(Xi) = sum_kl c_kl1 Xi_1^k Xi_2^l = X_x
	// T_y(Xi) = sum_kl c_kl2 Xi_1^k Xi_2^l = X_y
	// Newton's method 
	// find dXi such that
	// Xi + J_T(Xi)dXi = X

	auto eq = [epsilon](const double a, const double b) -> bool { return (fabs(a - b) <= epsilon * std::max(1.0, std::max(a, b))); };

	const Vector X({Xworld.x, Xworld.y});

	const auto f = [this, &X](const Vector& _Xi) -> std::pair<Scalar,Vector>
	{
		const Vector _T_Xi = getXworld(_Xi[0], _Xi[1]);
		const Matrix _J_Xi = getJacobian(_Xi[0], _Xi[1]);
		const Vector _diff =_T_Xi - X;

		return std::make_pair(norm(_diff), transpose(_J_Xi)*_diff / norm(_diff));
	};

	Vector Xi({0.0,0.0});

	for (size_t i=0;i<100;++i)
	{
		const Vector T_Xi = getXworld(Xi[0], Xi[1]);
		const Matrix inv_J = getInvJacobian(Xi[0], Xi[1]);

		if (eq(norm(T_Xi - X), 0.0)) { return std::make_pair(true, NodeRef(Xi[0], Xi[1])); }

		const Vector dXi = inv_J*(X - T_Xi);

		if (eq(norm(dXi), 0.0)) { return std::make_pair(false, NodeRef(Xi[0], Xi[1])); }

		Scalar alpha = std::get<1>( linesearchBisection(f, Xi, dXi, 1000, epsilon) );
		auto Xi_new = (Xi + alpha*dXi);

		for (size_t d=0;d<Xi.getSize();++d)
		{
			if (Xi_new[d] < -1.0) { alpha = Scalar(-1 + epsilon - Xi[d]) / Scalar(dXi[d]); }
			if (Xi_new[d] > 1.0) { alpha = Scalar(1 - epsilon - Xi[d]) / Scalar(dXi[d]);  }
		}

		if (eq(alpha, 0.0)) { return std::make_pair(false, NodeRef(Xi[0], Xi[1])); }

		Xi += alpha*dXi;

	}

	return std::make_pair(false, NodeRef(Xi[0], Xi[1]));
}
