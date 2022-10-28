/*
 * Analysis.hpp
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

#ifndef ANALYSIS_INCLUDE_HPP
#define ANALYSIS_INCLUDE_HPP

#include <LightFEM/Analysis/TrialFunction.hpp>
#include <LightFEM/Analysis/TestFunction.hpp>
#include <LightFEM/Analysis/VectorTrialFunction.hpp>
#include <LightFEM/Analysis/VectorTestFunction.hpp>
#include <LightFEM/Analysis/operators.hpp>

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp>
#include <LightFEM/Analysis/Measure/Measure.hpp>
#include <LightFEM/Analysis/Measure/GaussLobattoQuadrature.hpp>
#include <LightFEM/Analysis/Measure/SamplingFunction.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>
#include <LightFEM/Analysis/FunctionSpace/PnFunctionSpace.hpp>
#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>

#include <LightFEM/Analysis/Forms/LinearForm.hpp>
#include <LightFEM/Analysis/Forms/VectorLinearForm.hpp>
#include <LightFEM/Analysis/Forms/BilinearForm.hpp>
#include <LightFEM/Analysis/Forms/DiagonalBilinearForm.hpp>
#include <LightFEM/Analysis/Forms/VectorBilinearForm.hpp>

#endif // ANALYSIS_INCLUDE_HPP
