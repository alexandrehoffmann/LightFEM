/*
 * LightFEM-MPI.hpp
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

#ifndef LIGHT_FEM_MPI_HPP
#define LIGHT_FEM_MPI_HPP

#include <mpi.h>

#include <LightFEM/Analysis/Forms/MatrixEntry.hpp>

#include <LightFEM/Expression/ElementWise/ElementWiseFunction.hpp>
#include <LightFEM/Expression/Function/Function.hpp>

namespace LightFEM_MPI
{

void init();
void finalize();

void allReduce(const ScalarField& local_f, ScalarField& f, MPI_Comm comm);
void allReduce(const CpxScalarField& local_f, CpxScalarField& f, MPI_Comm comm);

////////////////////////////////////////////////////////////////////////

extern MPI_Datatype MPI_MATRIX_ENTRY;
extern MPI_Datatype MPI_CPX_MATRIX_ENTRY;
extern MPI_Datatype MPI_ELEMENT_WISE_SCALAR_FIELD;
extern MPI_Datatype MPI_CPX_ELEMENT_WISE_SCALAR_FIELD;

extern MPI_Op MPI_MATRIX_ENTRY_SUM;
extern MPI_Op MPI_CPX_MATRIX_ENTRY_SUM;
extern MPI_Op MPI_ELEMENT_WISE_SCALAR_FIELD_SUM;
extern MPI_Op MPI_CPX_ELEMENT_WISE_SCALAR_FIELD_SUM;

////////////////////////////////////////////////////////////////////////

void sum_MPI_MatrixEntry(void* in, void* inout, int* len, MPI_Datatype* MPI_matrixEntry);

void defineMatrixEntry();

////////////////////////////////////////////////////////////////////////

void sum_MPI_CpxMatrixEntry(void* in, void* inout, int* len, MPI_Datatype* MPI_matrixEntry);

void defineCpxMatrixEntry();

////////////////////////////////////////////////////////////////////////

void sum_MPI_ElementWiseScalarField(void* in, void* inout, int* len, MPI_Datatype* MPI_elementWiseScalarField);

void defineElementWiseScalarField();

////////////////////////////////////////////////////////////////////////

void sum_MPI_CpxElementWiseScalarField(void* in, void* inout, int* len, MPI_Datatype* MPI_elementWiseScalarField);

void defineCpxElementWiseScalarField();

////////////////////////////////////////////////////////////////////////

}

#endif // LIGHT_FEM_MPI_HPP
