/*
 * LightFEM-MPI.cpp
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

#include <LightFEM/LightFEM_MPI.hpp>

namespace LightFEM_MPI
{

MPI_Datatype MPI_MATRIX_ENTRY;
MPI_Datatype MPI_CPX_MATRIX_ENTRY;
MPI_Datatype MPI_ELEMENT_WISE_SCALAR_FIELD;
MPI_Datatype MPI_CPX_ELEMENT_WISE_SCALAR_FIELD;

MPI_Op MPI_MATRIX_ENTRY_SUM;
MPI_Op MPI_CPX_MATRIX_ENTRY_SUM;
MPI_Op MPI_ELEMENT_WISE_SCALAR_FIELD_SUM;
MPI_Op MPI_CPX_ELEMENT_WISE_SCALAR_FIELD_SUM;

void init()
{
	defineMatrixEntry();
	defineCpxMatrixEntry();
	
	defineElementWiseScalarField();
	defineCpxElementWiseScalarField();
	
	MPI_Op_create(static_cast<MPI_User_function*>(sum_MPI_MatrixEntry), 1, &MPI_MATRIX_ENTRY_SUM);
	MPI_Op_create(static_cast<MPI_User_function*>(sum_MPI_CpxMatrixEntry), 1, &MPI_CPX_MATRIX_ENTRY_SUM);
	
	MPI_Op_create(static_cast<MPI_User_function*>(sum_MPI_ElementWiseScalarField), 1, &MPI_ELEMENT_WISE_SCALAR_FIELD_SUM);
	MPI_Op_create(static_cast<MPI_User_function*>(sum_MPI_CpxElementWiseScalarField), 1, &MPI_CPX_ELEMENT_WISE_SCALAR_FIELD_SUM);
}

void finalize()
{
	MPI_Type_free(&MPI_MATRIX_ENTRY);
	MPI_Type_free(&MPI_CPX_MATRIX_ENTRY);
	MPI_Type_free(&MPI_ELEMENT_WISE_SCALAR_FIELD);
	MPI_Type_free(&MPI_CPX_ELEMENT_WISE_SCALAR_FIELD);
	
	
	MPI_Op_free(&MPI_MATRIX_ENTRY_SUM);
	MPI_Op_free(&MPI_CPX_MATRIX_ENTRY_SUM);
	MPI_Op_free(&MPI_ELEMENT_WISE_SCALAR_FIELD_SUM);
	MPI_Op_free(&MPI_CPX_ELEMENT_WISE_SCALAR_FIELD_SUM);
}

void allReduce(const ScalarField& local_f, ScalarField& f, MPI_Comm comm)
{
	if (local_f.getMesh() != f.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	
	MPI_Allreduce(local_f.data(), f.data(), f.getMesh()->getNElem(), MPI_ELEMENT_WISE_SCALAR_FIELD, MPI_ELEMENT_WISE_SCALAR_FIELD_SUM, comm);
	
	for (size_t e=0;e<f.getMesh()->getNElem();++e) 
	{
		f[e].setElement(local_f[e].getElement());
	}
}

void allReduce(const CpxScalarField& local_f, CpxScalarField& f, MPI_Comm comm)
{
	if (local_f.getMesh() != f.getMesh()) { throw std::invalid_argument("Both functions must be defined on the same mesh."); }
	
	MPI_Allreduce(local_f.data(), f.data(), f.getMesh()->getNElem(), MPI_CPX_ELEMENT_WISE_SCALAR_FIELD, MPI_CPX_ELEMENT_WISE_SCALAR_FIELD_SUM, comm);
	
	for (size_t e=0;e<f.getMesh()->getNElem();++e) 
	{
		f[e].setElement(local_f[e].getElement());
	}
}

////////////////////////////////////////////////////////////////////////

void sum_MPI_MatrixEntry(void *in, void *inout, int *len, MPI_Datatype*)
{
	MatrixEntry* invals    = static_cast<MatrixEntry*>(in);
	MatrixEntry* inoutvals = static_cast<MatrixEntry*>(inout);

	for (int i=0; i<*len; i++) 
	{
		inoutvals[i].row   += invals[i].row;
		inoutvals[i].col   += invals[i].col;
		inoutvals[i].value += invals[i].value;
	}
}

void defineMatrixEntry() 
{
	int          blocklens[3] = {1, 1, 1};
	MPI_Datatype types[3]     = {MPI_INT, MPI_INT, MPI_DOUBLE};
	MPI_Aint     disps[3]     = { offsetof(MatrixEntry, row),  offsetof(MatrixEntry, col), offsetof(MatrixEntry, value)};
	
	for (size_t i=1;i<3;++i) { disps[i] -= disps[0]; }
	disps[0] = 0;
	
	MPI_Type_create_struct(3, blocklens, disps, types, &MPI_MATRIX_ENTRY);
	MPI_Type_commit( &MPI_MATRIX_ENTRY );
}

////////////////////////////////////////////////////////////////////////

void sum_MPI_CpxMatrixEntry(void *in, void *inout, int *len, MPI_Datatype*)
{
	CpxMatrixEntry* invals    = static_cast<CpxMatrixEntry*>(in);
	CpxMatrixEntry* inoutvals = static_cast<CpxMatrixEntry*>(inout);

	for (int i=0; i<*len; i++) 
	{
		inoutvals[i].row   += invals[i].row;
		inoutvals[i].col   += invals[i].col;
		inoutvals[i].value += invals[i].value;
	}
}

void defineCpxMatrixEntry() 
{
	int          blocklens[3] = {1, 1, 1};
	MPI_Datatype types[3]     = {MPI_INT, MPI_INT, MPI_CXX_DOUBLE_COMPLEX};
	MPI_Aint     disps[3]     = { offsetof(CpxMatrixEntry, row),  offsetof(CpxMatrixEntry, col), offsetof(CpxMatrixEntry, value)};
	
	for (size_t i=1;i<3;++i) { disps[i] -= disps[0]; }
	disps[0] = 0;
	
	MPI_Type_create_struct(3, blocklens, disps, types, &MPI_CPX_MATRIX_ENTRY);
	MPI_Type_commit( &MPI_CPX_MATRIX_ENTRY );
}

////////////////////////////////////////////////////////////////////////

void sum_MPI_ElementWiseScalarField(void* in, void* inout, int* len, MPI_Datatype*)
{
	ElementWiseScalarField* invals    = static_cast<ElementWiseScalarField*>(in);
	ElementWiseScalarField* inoutvals = static_cast<ElementWiseScalarField*>(inout);
	
	for (int i=0; i<*len; i++) 
	{
		inoutvals[i].setElement(invals[i].getElement());
		inoutvals[i] += invals[i];
	}
}

void defineElementWiseScalarField()
{
	MPI_Type_contiguous(sizeof(ElementWiseScalarField), MPI_BYTE, &MPI_ELEMENT_WISE_SCALAR_FIELD);
	MPI_Type_commit(&MPI_ELEMENT_WISE_SCALAR_FIELD);
}

////////////////////////////////////////////////////////////////////////

void sum_MPI_CpxElementWiseScalarField(void* in, void* inout, int* len, MPI_Datatype*)
{
	CpxElementWiseScalarField* invals    = static_cast<CpxElementWiseScalarField*>(in);
	CpxElementWiseScalarField* inoutvals = static_cast<CpxElementWiseScalarField*>(inout);
	
	for (int i=0; i<*len; i++) 
	{
		inoutvals[i].setElement(invals[i].getElement());
		inoutvals[i] += invals[i];
	}
}

void defineCpxElementWiseScalarField()
{
	MPI_Type_contiguous(sizeof(CpxElementWiseScalarField), MPI_BYTE, &MPI_CPX_ELEMENT_WISE_SCALAR_FIELD);
	MPI_Type_commit(&MPI_CPX_ELEMENT_WISE_SCALAR_FIELD);
}

////////////////////////////////////////////////////////////////////////

}
