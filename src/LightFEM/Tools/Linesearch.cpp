/*
 * Linesearch.cpp
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

#include <LightFEM/Tools/Linesearch.hpp>

#include <LightFEM/Expression/LinAlg/operators.hpp>
//#include <iostream>

std::pair< bool, Scalar > linesearchBisection(const std::function< std::pair<Scalar, Vector> (const Vector& x)>& f, const Vector& x0, const Vector& d, const size_t maxIt, const double epsilon)
{
	auto eq = [epsilon](const double a, const double b) -> bool { return (fabs(a - b) <= epsilon * std::max(1.0, std::max(a, b))); };
	auto leq = [&eq](const double a, const double b) -> bool { return (a < b) or eq(a,b); };
	//auto s_eq = [&eq](const Scalar a, const Scalar b) -> bool { return eq(a.eval(), b.eval()); };
	//auto s_leq = [&leq](const Scalar a, const Scalar b) -> bool { return leq(a.eval(), b.eval()); };


	const Scalar c1(1.0e-4);
	const Scalar c2(0.9);

	Scalar alpha_max(std::numeric_limits< double >::infinity());
	Scalar alpha_min(0.0);
	Scalar alpha(1.0);

	const auto [fx0, grad_fx0] = f(x0);

	const Scalar m0 = inner(d, grad_fx0);

	for (size_t it=0;it<maxIt;++it)
	{
		if (std::isfinite(alpha_max) and eq(alpha_min, alpha_max)) { break; }

		const Vector xk = x0 + alpha*d;
		const auto [fxk, grad_fxk] = f(xk);
		const Scalar mk = inner(d, grad_fxk);


		const bool isArmijoSatisfied = leq(fxk, fx0 + c1*alpha*m0);
		const bool isCurvatureSatisfied = leq(-mk, -c2*m0);

		if (isArmijoSatisfied and isCurvatureSatisfied)
		{
			return std::make_pair(true, alpha);
		}
		if (not isArmijoSatisfied) // step too long
		{
			alpha_max = alpha;
			alpha = Scalar(0.5)*(alpha_max + alpha_min);
		}
		else if (not std::isfinite(alpha_max) and not isCurvatureSatisfied) // step too short
		{
			alpha *= Scalar(2.0);
		}
		else if (not isCurvatureSatisfied) // step too short
		{
			alpha_min = alpha;
			alpha = Scalar(0.5)*(alpha_max + alpha_min);
		}
	}
	return std::make_pair(false, alpha);
}
