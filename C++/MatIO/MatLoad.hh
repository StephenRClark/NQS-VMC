// ------------------------------------------------------------------------------------------------------------
// MatLoad.hh - contains the processes necessary to convert from loaded .mat files containing NQS Modifiers 
//					to the C++ equivalents and back. Necessary for interfacing with the Matlab code.
// ------------------------------------------------------------------------------------------------------------

#ifndef MATLOAD_HH
#define MATLOAD_HH

#include <iostream>
#include <string.h>
#include <Eigen/Dense>
#include "mat.h"
#include "matrix.h"
#include "Hilbert/Hilbert.hh"
#include "Graph/Graph.hh"
#include "Operator/Operator.hh"
#include "Hamiltonian/Hamiltonian.hh"
#include "Sampler/MCSampler.hh"
#include "Optimiser/StochReconfig.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// ListField will list the fields available in the mxArray pointed to by StrPtr.
	void ListField(const mxArray* StrPtr);
	// mex2eig_vec is a basic conversion of a mxArray vector to an Eigen vector.
	VectorXd mex2eig_vec(const mxArray* matvec);
	// ParamExtract will search the file for the given variable name and attempt to output an Eigen vector.
	VectorXd ParamExtract(const char* filename, const char* paramname);
	// HilbertExtract will search the file for the given struct name and attempt to construct a Hilbert object.
	Hilbert* HilbertExtract(const char* filename, const char* structname);
	// GraphExtract will search the file for the given struct name and attempt to construct a Graph object.
	Graph* GraphExtract(const char* filename, const char* structname);
	// OperatorExtract will search the file for the given struct name and attempt to construct an Operator object.
	Operator* OperatorExtract(const char* filename, const char* structname);
	// HamiltionianExtract will file for the given struct name and attempt to construct a Hamiltonian object.
	Hamiltonian* HamiltonianExtract(const char* filename, const char* structname);
	// ReferenceExtract will search the file for the given struct name and attempt to construct a Reference object.
	// - Requires a Hilbert instance to perform.
	Reference* ReferenceExtract(const char* filename, const char* structname, Hilbert* hlb);
	// ModifierExtract will search the file for the given struct name and attempt to construct a Modifier object.
	// - Requires a Hilbert instance to perform.
	Modifier* ModifierExtract(const char* filename, const char* structname, Hilbert* hlb);
	// AnsatzExtract will search the file for the given struct name and attempt to construct an Ansatz object.
	Ansatz* AnsatzExtract(const char* filename, const char* structname);
	// SamplerExtract will search the file for the given struct name and attempt to construct a MCSampler object.
	// - Requires a Hilbert instance to perform.
	MCSampler* SamplerExtract(const char* filename, const char* structname, Hilbert* hlb);
	// SRExtract will search the file for the given struct name and attempt to construct a StochReconfig object.
	StochReconfig* SRExtract(const char* filename, const char* structname);

	// VectorMatSave will take a vector (Eig or std) and save/append it to the given filename with the given variable name.
	// - Eigen vector input version:
	void VectorMatSave(VectorXd EigVec, const char* filename, const char* varname);
	// - std vector input version:
	void VectorMatSave(vector<double> StdVec, const char* filename, const char* varname);
	// MatrixMatSave will take a matrix (Eig or std) and save/append it to the given filename with the given variable name.
	// - Eigen matrix input version:
	void MatrixMatSave(MatrixXd EigMat, const char* filename, const char* varname);
	// - std vector input version:
	void MatrixMatSave(vector<vector<double>> StdMat, const char* filename, const char* varname);
	// ScalarMatSave will take a double input and save/append it to the given filename with the given variable name.
	void ScalarMatSave(double StdDbl, const char* filename, const char* varname);
}

#endif