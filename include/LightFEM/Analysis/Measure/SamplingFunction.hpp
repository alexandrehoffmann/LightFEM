/*
 * SamplingFunction.hpp
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

#ifndef SAMPLING_FUNCTION_HPP
#define SAMPLING_FUNCTION_HPP

#include <LightFEM/Analysis/Measure/Measure.hpp>

class SamplingFunction : public Measure
{
public:
	SamplingFunction(const Mesh* mesh, const NodeWorld& X);
	SamplingFunction(const Mesh* mesh, const NodeWorld& X, double weight);
	SamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X);
	SamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< double >& weights);
private:
	void init(const Mesh* mesh, const std::vector<NodeWorld>& X, const std::vector< double >& weights);
public:
	static double epsilon;
};

////////////////////////////////////////////////////////////////////////

class CpxSamplingFunction : public CpxMeasure
{
public:
	CpxSamplingFunction(const Mesh* mesh, const NodeWorld& X);
	CpxSamplingFunction(const Mesh* mesh, const NodeWorld& X, const std::complex< double > weight);
	CpxSamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X);
	CpxSamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< std::complex< double > >& weights);
private:
	void init(const Mesh* mesh, const std::vector<NodeWorld>& X, const std::vector< std::complex< double > >& weights);
public:
	static double epsilon;
};

#undef UNARY_TEMPLATE_MACRO

#endif // SAMPLING_FUNCTION_HPP
