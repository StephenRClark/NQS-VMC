 // ------------------------------------------------------------------------------------------------------------
// MatLoad.cpp - contains the processes necessary to convert from loaded .mat files containing NQS Modifiers 
//					to the C++ equivalents and back. Necessary for interfacing with the Matlab code.
// ------------------------------------------------------------------------------------------------------------

#ifndef MATLOAD_CPP
#define MATLOAD_CPP

#include "MatLoad.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// ListField will list the fields available in the mxArray pointed to by StrPtr.
	void ListField(const mxArray* StrPtr)
	{
		int nfield = mxGetNumberOfFields(StrPtr);
		const char* classname;
		classname = mxGetClassName(StrPtr);
		std::cout << "Input mxArray is of class " << classname << " and has " << nfield << " fields." << endl;
		for (int f = 0; f < nfield; f++)
		{
			const char* fieldname;
			fieldname = mxGetFieldNameByNumber(StrPtr, f);
			std::cout << "Field " << f << ": " << fieldname << endl;
		}
		return;
	}
	// FindField is a function that will look for the field in the provided mxArray with the provided fieldname.
	const mxArray* FindField(const mxArray* StrPtr, const char* fieldname)
	{
		size_t numfield = mxGetNumberOfElements(StrPtr);
		const mxArray* FldPtr = NULL;
		for (int f = 0; f < numfield; f++)
		{
			FldPtr = mxGetField(StrPtr, f, fieldname);
			if (FldPtr != NULL)
			{
				break;
			}
		}
		if (FldPtr == NULL)
		{
			cout << "Could not find field " << fieldname << " - returning null pointer." << endl;
		}
		return FldPtr;
	}
	// mex2eig_vec is a basic conversion of a mxArray vector to an Eigen vector.
	VectorXd mex2eig_vec(const mxArray* matvec)
	{		
		double* pr;
		pr = mxGetPr(matvec);
		int size = (int)mxGetNumberOfElements(matvec);
		VectorXd eigvec = VectorXd::Zero(size);
		for (int p = 0; p < size; p++)
		{
			eigvec(p) = pr[p];
		}
		return eigvec;
	}
	// mex2std_dvec is a basic conversion of a mxArray vector to an std double vector.
	vector<double> mex2std_dvec(const mxArray* matvec)
	{
		double* pr;
		pr = mxGetPr(matvec);
		int size = (int)mxGetNumberOfElements(matvec);
		vector<double> stdvec(size);
		for (int p = 0; p < size; p++)
		{
			stdvec[p] = pr[p];
		}
		return stdvec;
	}
	// mex2std_ivec is a basic conversion of a mxArray vector to an std int vector.
	vector<int> mex2std_ivec(const mxArray* matvec)
	{
		double* pr;
		pr = mxGetPr(matvec);
		int size = (int)mxGetNumberOfElements(matvec);
		vector<int> stdvec(size);
		for (int p = 0; p < size; p++)
		{
			stdvec[p] = (int)pr[p];
		}
		return stdvec;
	}
	// mex2dbl is a basic conversion of a mxArray value to a double value.
	double mex2dbl(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		if (size > 1)
		{
			std::cout << "Input mxArray has more than one element - output will only contain first element." << endl;
		}
		double* pr;
		pr = mxGetPr(matvec);
		double val = pr[0];
		return val;
	}
	// mex2int is a basic conversion of a mxArray value to an int value.
	int mex2int(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		if (size > 1)
		{
			std::cout << "Input mxArray has more than one element - output will only contain first element." << endl;
		}
		double* pr;
		pr = mxGetPr(matvec);		
		int val = (int)pr[0];
		return val;
	}
	// mex2st_dmat is a conversion of a mxArray matrix to a nested double vector matrix. 
	// - Column major variant (columns of Matlab matrix as individual vectors).
	vector<vector<double>> mex2std_dmat_c(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		const size_t* dimensions = mxGetDimensions(matvec);
		vector<vector<double>> mat(dimensions[1]);
		double* pr;
		pr = mxGetPr(matvec);
		for (size_t m = 0; m < dimensions[1]; m++)
		{
			vector<double> col(dimensions[0]);
			for (size_t n = 0; n < dimensions[0]; n++)
			{
				col[n] = pr[n + (m * dimensions[0])];
			}
			mat[m] = col;
		}
		return mat;
	}
	// mex2st_dmat is a conversion of a mxArray matrix to a nested double vector matrix. 
	// - Row major variant (rows of Matlab matrix as individual vectors).
	vector<vector<double>> mex2std_dmat_r(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		const size_t* dimensions = mxGetDimensions(matvec);
		vector<vector<double>> mat(dimensions[0]);
		double* pr;
		pr = mxGetPr(matvec);
		for (size_t m = 0; m < dimensions[0]; m++)
		{
			vector<double> col(dimensions[1]);
			for (size_t n = 0; n < dimensions[1]; n++)
			{
				col[n] = pr[m + (n * dimensions[1])];
			}
			mat[m] = col;
		}
		return mat;
	}
	// mex2st_imat is a conversion of a mxArray matrix to a nested int vector matrix. 
	// - Column major variant (columns of Matlab matrix as individual vectors).
	vector<vector<int>> mex2std_imat_c(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		const size_t* dimensions = mxGetDimensions(matvec);
		vector<vector<int>> mat(dimensions[1]);
		double* pr;
		pr = mxGetPr(matvec);
		for (size_t m = 0; m < dimensions[1]; m++)
		{
			vector<int> col(dimensions[0]);
			for (size_t n = 0; n < dimensions[0]; n++)
			{
				col[n] = (int)pr[n + (m * dimensions[0])];
			}
			mat[m] = col;
		}		
		return mat;
	}
	// mex2std_imat is a conversion of a mxArray matrix to a nested int vector matrix. 
	// - Row major variant (rows of Matlab matrix as individual vectors).
	vector<vector<int>> mex2std_imat_r(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		const size_t* dimensions = mxGetDimensions(matvec);
		vector<vector<int>> mat(dimensions[0]);
		double* pr;
		pr = mxGetPr(matvec);
		for (size_t m = 0; m < dimensions[0]; m++)
		{
			vector<int> col(dimensions[1]);
			for (size_t n = 0; n < dimensions[1]; n++)
			{
				col[n] = (int)pr[m + (n * dimensions[1])];
			}
			mat[m] = col;
		}
		return mat;
	}
	// mex2eig_mat is a conversion of a mxArray matrix to an Eigen matrix. 
	MatrixXd mex2eig_mat(const mxArray* matvec)
	{
		int size = (int)mxGetNumberOfElements(matvec);
		const size_t* dimensions = mxGetDimensions(matvec);
		MatrixXd mat = MatrixXd::Zero(dimensions[0],dimensions[1]);
		double* pr;
		pr = mxGetPr(matvec);
		for (size_t m = 0; m < dimensions[0]; m++)
		{
			for (size_t n = 0; n < dimensions[1]; n++)
			{
				mat(m,n) = pr[m + (n * dimensions[0])];
			}
		}
		return mat;
	}

	// ParamExtract will grab any loaded Params vector in the mat file and convert to an Eigen vector.
	VectorXd ParamExtract(const char* filename, const char* paramname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << " - returning null pointer. Check name is correct." << endl;
			std::abort();
		}
		const mxArray* MatParams;
		MatParams = matGetVariable(matfile, paramname);
		if (MatParams == NULL)
		{
			cerr << "Cannot find Params in .mat file." << endl;
			std::abort();
		}
		VectorXd EigParams = mex2eig_vec(MatParams);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		return EigParams;
	}

	// BoseConstruct will take a Bose Hilbert from Matlab (identified in HilbertConstruct) and give the C++ equivalent.
	Hilbert* BoseConstruct(const mxArray* hilbertstruct)
	{
		// Extract information common to all Bose subclasses.
		const char* fieldname;
		const mxArray* fieldptr;		
		Hilbert* h_;
		int N = 1;
		double* Nmat;
		fieldptr = FindField(hilbertstruct, "N");
		if (mxIsDouble(fieldptr))
		{
			Nmat = mxGetPr(fieldptr);
			N = (int)Nmat[0];
		}
		int Nb = 1;
		double* Nbmat;
		fieldptr = FindField(hilbertstruct, "Nb");
		if (mxIsDouble(fieldptr))
		{
			Nbmat = mxGetPr(fieldptr);
			Nb = (int)Nbmat[0];
		}
		int Nmax = 1;
		double* Nmaxmat;
		fieldptr = FindField(hilbertstruct, "Nmax");
		if (mxIsDouble(fieldptr))
		{
			Nmaxmat = mxGetPr(fieldptr);
			Nmax = (int)Nmaxmat[0];
		}		
		// Establish whether this is the fixed or variable total number case.		
		fieldname = "Sector";
		fieldptr = FindField(hilbertstruct, fieldname);
		if (mxIsEmpty(fieldptr)) // Empty corresponds to no Sector value - variable number.
		{
			h_= new Hilbert(N, Nb, Nmax, false); // Gives a Hilbert containing a BoseVN subclass pointer.
		}
		else
		{
			Nb = (int)mxGetPr(fieldptr)[0]; // Sector sets the number of bosons if not empty.
			h_ = new Hilbert(N, Nb, Nmax, true); // Gives a Hilbert containing a BoseFN subclass pointer.
		}
		return h_;
	}
	// FermConstruct will take a Ferm Hilbert from Matlab (identified in HilbertConstruct) and give the C++ equivalent.
	Hilbert* FermConstruct(const mxArray* hilbertstruct)
	{
		// Extract information common to all Bose subclasses.
		const char* fieldname;
		const mxArray* fieldptr;		
		int N = 1;
		double* Nmat;
		fieldptr = FindField(hilbertstruct, "N");
		if (mxIsDouble(fieldptr))
		{
			Nmat = mxGetPr(fieldptr);
			N = (int)Nmat[0];
		}
		int Nf = 1;
		double* Nfmat;
		fieldptr = FindField(hilbertstruct, "Nf");
		if (mxIsDouble(fieldptr))
		{
			Nfmat = mxGetPr(fieldptr);
			Nf = (int)Nfmat[0];
		}
		int Sz = 1;
		double* Szmat;
		fieldptr = FindField(hilbertstruct, "Sz");
		if (mxIsDouble(fieldptr))
		{
			Szmat = mxGetPr(fieldptr);
			Sz = (int)Szmat[0];
		}
		// Establish whether this is the fixed or variable total number case.		
		fieldname = "Sector";
		fieldptr = FindField(hilbertstruct, fieldname);
		double* Sector;
		char constraint;
		if (mxIsDouble(fieldptr))
		{
			Sector = mxGetPr(fieldptr);
			if (Sector[0] == 0) // No fixed particle number.
			{
				constraint = 'x';
			}
			else if (Sector[1] == 0) // No fixed spin projection.
			{
				constraint = 'n';
			}
			else // Fixed number and spin projection.
			{
				constraint = 's';
			}			
		}
		else // If no Sector information, assume unconstrained.
		{
			constraint = 'x';
		}
		Hilbert* h_ = new Hilbert(N, Nf, Sz, constraint);
		return h_;
	}
	// SpinConstruct will take a Spin Hilbert from Matlab (identified in HilbertConstruct) and give the C++ equivalent.
	Hilbert* SpinConstruct(const mxArray* hilbertstruct)
	{
		// Extract information common to all Bose subclasses.
		const char* fieldname;
		const mxArray* fieldptr;
		Hilbert* h_;
		int N = 1;
		double* Nmat;
		fieldptr = FindField(hilbertstruct, "N");
		if (mxIsDouble(fieldptr))
		{
			Nmat = mxGetPr(fieldptr);
			N = (int)Nmat[0];
		}
		double S = 0.5;
		double* Smat;
		fieldptr = FindField(hilbertstruct, "S");
		if (mxIsDouble(fieldptr))
		{
			Smat = mxGetPr(fieldptr);
			S = Smat[0];
		}
		int SzTotal;
		// Establish whether this is the fixed or variable total number case.		
		fieldname = "Sector";
		fieldptr = FindField(hilbertstruct, fieldname);
		if (mxIsEmpty(fieldptr)) // Empty corresponds to no Sector value - variable number.
		{
			SzTotal = 0; 
			h_ = new Hilbert(N, SzTotal, S, false); // Gives a Hilbert containing a SpinVZ subclass pointer.
		}
		else
		{
			double* Sector = mxGetPr(fieldptr);
			SzTotal = (int)Sector[0];
			h_ = new Hilbert(N, SzTotal, S, true); // Gives a Hilbert containing a SpinFZ subclass pointer.
		}
		return h_;
	}
	// HilbertConstruct takes a Hilbert Matlab struct and outputs a C++ equivalent.
	Hilbert* HilbertConstruct(const mxArray* structptr)
	{
		Hilbert* hlb = NULL;
		const mxArray* matvar;
		matvar = FindField(structptr, "Nb"); // Nb unique to Bose Hilbert objects.
		if (matvar != NULL)
		{
			hlb = BoseConstruct(structptr);
			return hlb;
		}
		matvar = FindField(structptr, "S"); // S unique to Spin Hilbert objects.
		if (matvar != NULL)
		{
			hlb = SpinConstruct(structptr);
			return hlb;
		}
		matvar = FindField(structptr, "Nf"); // Nf unique to Ferm Hilbert objects.
		if (matvar != NULL)
		{
			hlb = FermConstruct(structptr);
			return hlb;
		}
		return hlb;
	}
	// HilbertExtract will search the file for the given struct name and attempt to construct a Hilbert object.
	Hilbert* HilbertExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << " - returning null pointer. Check name is correct." << endl;
			std::abort();
		}
		// Each Hilbert type has a property unique to them - use this to work out which type.
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Hilbert* hlb = HilbertConstruct(propstruct);		
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (hlb == NULL)
		{
			cerr << "Unable to construct appropriate Hilbert object - check struct fields match." << endl;
			std::abort();
		}
		return hlb;
	}

	// GraphConstruct takes a Graph Matlab struct and outputs the C++ equivalent.
	Graph* GraphConstruct(const mxArray* structptr)
	{		
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Dim";
		matvar = FindField(structptr, fieldname);
		vector<int> dim = mex2std_ivec(matvar);		
		fieldname = "LVecs";
		matvar = FindField(structptr, fieldname);
		vector<vector<int>> lvecs = mex2std_imat_r(matvar);
		fieldname = "Bound";
		matvar = FindField(structptr, fieldname);
		vector<int> bound = mex2std_ivec(matvar);
		fieldname = "Bonds";
		matvar = FindField(structptr, fieldname);
		vector<vector<int>> bonds = mex2std_imat_c(matvar);
		for (int b = 0; b < bonds.size(); b++) // Need to convert from 1-based to 0-based.
		{
			for (int n = 0; n < bonds[b].size(); n++)
			{
				bonds[b][n] -= 1;
			}
		}
		fieldname = "SLInds";
		matvar = FindField(structptr, fieldname);
		vector<int> slinds = mex2std_ivec(matvar);
		fieldname = "Ntr";
		matvar = FindField(structptr, fieldname);
		int Ntr = mex2int(matvar);
		int nvecs = (int)lvecs.size();
		bool bondmapflag = true;
		if (Ntr == nvecs) // No need to use BondMapGen - this is only the case if the original did not use BondMapGen.
		{
			bondmapflag = false;
		}
		Graph* grp = new Graph(dim, lvecs, bound, bonds, slinds, bondmapflag);		
		if (grp == NULL)
		{
			cerr << "Unable to construct appropriate Graph object - check struct fields match." << endl;
			std::abort();
		}
		return grp;
	}
	// GraphExtract will search the file for the given struct name and attempt to construct a Graph object.
	Graph* GraphExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Graph* grp = GraphConstruct(propstruct);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (grp == NULL)
		{
			cerr << "Unable to construct appropriate Graph object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Graph constructed successfully." << endl;
		return grp;
	}

	// GetOpName will take a mxChar* pointer for a Matlab character vector and give the corresponding Operator signature
	const char* GetOpName(const char* funcname)
	{
		vector<const char*> checkname; // May replace these lists with externally defined ones later.
		vector<const char*> returnname;
		const char* endname;
		endname = "blank";
		// Energy gradient operators.
		checkname.push_back("SR_EnLogDerivCorr"); returnname.push_back("EdP");
		checkname.push_back("SR_LogDerivCorr"); returnname.push_back("dPdQ");
		checkname.push_back("TE_EnLocSqCorr"); returnname.push_back("Ensq");
		// Bosonic operators - update as more are incorporated.
		checkname.push_back("BpBm_OpMatEls"); returnname.push_back("BpBm");		
		checkname.push_back("NiNj_Bose_CfgVal"); returnname.push_back("NiNj");
		checkname.push_back("VarN_Bose_CfgVal"); returnname.push_back("VarN");
		checkname.push_back("OccFrac_Bose_CfgVal"); returnname.push_back("OcFr");
		checkname.push_back("DbHl_Bose_CfgVal"); returnname.push_back("DbHl");
		for (int o = 0; o < checkname.size(); o++)
		{
			if (strcmp(checkname[o], funcname) == 0)
			{
				endname = returnname[o];
				break;
			}
		}
		return endname;
	}
	// OperatorConstruct takes an Operator Matlab struct and outputs the C++ equivalent.
	Operator* OperatorConstruct(const mxArray* structptr)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Hilbert";
		matvar = FindField(structptr, fieldname);
		Hilbert* hlb = HilbertConstruct(matvar);
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* grp = GraphConstruct(matvar);
		fieldname = "FuncHandle";
		matvar = FindField(structptr, fieldname);
		const char* cfunc;
		cfunc = "blank";
		if (mxIsChar(matvar))
		{
			const char* funcname = mxArrayToString(matvar);
			cfunc = GetOpName(funcname);
		}
		else
		{
			cerr << "Cannot find FuncHandle character vector - check fields of the input struct." << endl;
			std::abort();
		}
		Operator* op_ = new Operator(hlb, grp, cfunc);
		return op_;
	}
	// OperatorExtract will search the file for the given struct name and attempt to construct an Operator object.
	Operator* OperatorExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Operator* op_ = OperatorConstruct(propstruct);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (op_ == NULL)
		{
			cerr << "Unable to construct appropriate Operator object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Operator constructed successfully." << endl;
		return op_;
	}

	// HamiltonianConstruct takes a Hamiltonian Matlab struct and outputs the C++ equivalent.
	Hamiltonian* HamiltonianConstruct(const mxArray* structptr)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "HParams";
		matvar = FindField(structptr, fieldname);
		vector<double> hparams = mex2std_dvec(matvar);
		fieldname = "Operator";
		matvar = FindField(structptr, fieldname); // Operators will be in a cell array.
		size_t numop = mxGetNumberOfElements(matvar);
		const mxArray* cellvar;
		vector<Operator*> hamop_(numop);
		for (int o = 0; o < numop; o++)
		{			
			cellvar = mxGetCell(matvar, o);
			hamop_[o] = OperatorConstruct(cellvar);
		}
		Hamiltonian* hm_ = new Hamiltonian(hparams, hamop_);
		return hm_;
	}
	// HamiltonianExtract will search the file for the given struct name and attempt to construct a Hamiltonian object.
	Hamiltonian* HamiltonianExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Hamiltonian* hm_ = HamiltonianConstruct(propstruct);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (hm_ == NULL)
		{
			cerr << "Unable to construct appropriate Hamiltonian object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Hamiltonian constructed successfully." << endl;
		return hm_;
	}

	// ReferenceConstruct takes a Reference Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Reference* ReferenceConstruct(const mxArray* structptr, Hilbert* hlb)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Type";
		matvar = FindField(structptr, fieldname);
		const char* reftype = mxArrayToString(matvar);
		Reference* ref = NULL;
		if (strcmp(reftype, "BECR") == 0)
		{
			fieldname = "SPO";
			matvar = FindField(structptr, fieldname);
			VectorXd SPO = mex2eig_vec(matvar);
			ref = new BECR(hlb, SPO);
		}
		// Other cases to be added later.
		return ref;
	}
	// ReferenceExtract will search the file for the given struct name and attempt to construct a Reference object.
	// - Requires a Hilbert instance to perform.
	Reference* ReferenceExtract(const char* filename, const char* structname, Hilbert* hlb)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Reference* ref = ReferenceConstruct(propstruct, hlb);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (ref == NULL)
		{
			cerr << "Unable to construct appropriate Reference object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Reference constructed successfully." << endl;
		return ref;
	}

	// NQSSHConstruct takes a NQSSHTI Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* NQSSHConstruct(const mxArray* structptr, Hilbert* hlb)
	{ // Default case is translation invariant - C++ non-translation invariant class is NQSSHns.
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Nv";
		matvar = FindField(structptr, fieldname);
		int Nv = mex2int(matvar);
		fieldname = "Alpha";
		matvar = FindField(structptr, fieldname);
		int Alpha = mex2int(matvar);
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "HDim";
		matvar = FindField(structptr, fieldname);
		int HDim = mex2int(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		// Construct a new Hilbert, only really needed for initialisation.
		Hilbert* TempHilbert = new Hilbert(Nv, Nv, (HDim - 1), true);
		Modifier* nqs_ = new NQSSH(TempHilbert, GraphObj, Alpha, Params);		
		nqs_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		nqs_->SetParamCap(ParamCap);
		return nqs_;
	}
	// NQSNHConstruct takes a NQSNHTI Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* NQSNHConstruct(const mxArray* structptr, Hilbert* hlb)
	{ // Default case is translation invariant - C++ non-translation invariant class is NQSNHns.
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Nv";
		matvar = FindField(structptr, fieldname);
		int Nv = mex2int(matvar);
		fieldname = "Alpha";
		matvar = FindField(structptr, fieldname);
		int Alpha = mex2int(matvar);
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "HDim";
		matvar = FindField(structptr, fieldname);
		int HDim = mex2int(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		// Construct a new Hilbert, only really needed for initialisation.
		Hilbert* TempHilbert = new Hilbert(Nv, Nv, (HDim - 1), true);
		Modifier* nqs_ = new NQSNH(TempHilbert, GraphObj, Alpha, Params);		
		nqs_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		nqs_->SetParamCap(ParamCap);
		return nqs_;
	}
	// NQSMHConstruct takes a NQSMHTI Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* NQSMHConstruct(const mxArray* structptr, Hilbert* hlb)
	{ // Default case is translation invariant - C++ non-translation invariant class is NQSMHns.
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Nv";
		matvar = FindField(structptr, fieldname);
		int Nv = mex2int(matvar);
		fieldname = "Alpha";
		matvar = FindField(structptr, fieldname);
		int Alpha = mex2int(matvar);
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "HDim";
		matvar = FindField(structptr, fieldname);
		int HDim = mex2int(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		// Construct a new Hilbert, only really needed for initialisation.
		Hilbert* TempHilbert = new Hilbert(Nv, Nv, (HDim - 1), true);
		Modifier* nqs_ = new NQSMH(TempHilbert, GraphObj, Alpha, Params);		
		nqs_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		nqs_->SetParamCap(ParamCap);
		return nqs_;
	}
	// NQSOHConstruct takes a NQSOHTI Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* NQSOHConstruct(const mxArray* structptr, Hilbert* hlb)
	{ // Default case is translation invariant - C++ non-translation invariant class is NQSOHns.
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Alpha";
		matvar = FindField(structptr, fieldname);
		int Alpha = mex2int(matvar);
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "VDim";
		matvar = FindField(structptr, fieldname);
		int VDim = mex2int(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		Modifier* nqs_ = new NQSOH(hlb, GraphObj, Alpha, Params);		
		nqs_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		nqs_->SetParamCap(ParamCap);
		return nqs_;
	}
	// JastConstruct takes a Jast Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* JastConstruct(const mxArray* structptr, Hilbert* hlb)
	{ // Default case is translation invariant - C++ non-translation invariant class is Jastns.
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		Modifier* j_ = new Jast(hlb, GraphObj, Params);		
		j_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		j_->SetParamCap(ParamCap);
		return j_;
	}
	// NNMBConstruct takes a NNMB Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* NNMBConstruct(const mxArray* structptr, Hilbert* hlb)
	{ 
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Graph";
		matvar = FindField(structptr, fieldname);
		Graph* GraphObj = GraphConstruct(matvar);
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		Params = Params.array() * OptInds.array();
		Modifier* mbc_ = new NNMB(hlb, GraphObj, Params);		
		mbc_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		mbc_->SetParamCap(ParamCap);
		return mbc_;
	}
	// GutzConstruct takes a GutzB Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* GutzConstruct(const mxArray* structptr, Hilbert* hlb)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Params";
		matvar = FindField(structptr, fieldname);
		VectorXd Params = mex2eig_vec(matvar);
		Modifier* g_ = new Gutz(hlb, Params(0), true);
		fieldname = "OptInds";
		matvar = FindField(structptr, fieldname);
		VectorXd OptInds = mex2eig_vec(matvar);
		g_->OptIndLoad(OptInds);
		fieldname = "ParamCap";
		matvar = FindField(structptr, fieldname);
		double ParamCap = mex2dbl(matvar);
		g_->SetParamCap(ParamCap);
		return g_;
	}
	// ModifierConstruct takes a Modifier Matlab struct and a pre-made Hilbert and outputs the C++ equivalent.
	Modifier* ModifierConstruct(const mxArray* structptr, Hilbert* hlb)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Type";
		matvar = FindField(structptr, fieldname);
		const char* modtype = mxArrayToString(matvar);
		Modifier* mod = NULL;
		if (strcmp(modtype, "NQSSHTI") == 0)
		{
			mod = NQSSHConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "NQSNHTI") == 0)
		{
			mod = NQSNHConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "NQSMHTI") == 0)
		{
			mod = NQSMHConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "NQSOHTI") == 0)
		{
			mod = NQSOHConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "GutzB") == 0)
		{
			mod = GutzConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "Jast") == 0)
		{
			mod = JastConstruct(structptr, hlb);
		}
		else if (strcmp(modtype, "NNMB") == 0)
		{
			mod = NNMBConstruct(structptr, hlb);
		}
		// Other cases to be added later.
		return mod;
	}
	// ModifierExtract will search the file for the given struct name and attempt to construct a Modifier object.
	// - Requires a Hilbert instance to perform.
	Modifier* ModifierExtract(const char* filename, const char* structname, Hilbert* hlb)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Modifier* mod = ModifierConstruct(propstruct, hlb);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (mod == NULL)
		{
			cerr << "Unable to construct appropriate Modifier object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Modifier constructed successfully." << endl;
		return mod;
	}

	// AnsatzConstruct takes an Ansatz Matlab struct and outputs the C++ equivalent.
	Ansatz* AnsatzConstruct(const mxArray* structptr)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Hilbert";
		matvar = FindField(structptr, fieldname);
		Hilbert* hlb = HilbertConstruct(matvar); // Use generated Hilbert to initialise Reference.
		cout << "Ansatz Hilbert constructed successfully." << endl;
		fieldname = "Reference";
		matvar = FindField(structptr, fieldname); 
		// Reference and Hilbert must be compatible in Matlab code, so shouldn't run into issues using Hilbert here.
		Reference* ref = ReferenceConstruct(matvar, hlb);
		cout << "Ansatz Reference constructed successfully." << endl;
		fieldname = "Modifier";
		matvar = FindField(structptr, fieldname);
		size_t nummod = mxGetNumberOfElements(matvar);// Operators will be in a cell array.
		const mxArray* cellvar;
		vector<Modifier*> mod(nummod);
		for (int o = 0; o < nummod; o++)
		{
			cellvar = mxGetCell(matvar, o);
			mod[o] = ModifierConstruct(cellvar,hlb);
		}
		cout << "Ansatz Modifiers constructed successfully." << endl;
		Ansatz* ans = new Ansatz(ref,mod,hlb);
		return ans;
	}
	// AnsatzExtract will search the file for the given struct name and attempt to construct an Ansatz object.
	Ansatz* AnsatzExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		Ansatz* ans = AnsatzConstruct(propstruct);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (ans == NULL)
		{
			cerr << "Unable to construct appropriate Ansatz object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New Ansatz constructed successfully." << endl;
		return ans;
	}
	
	// SamplerConstruct takes a Sampler Matlab struct and outputs the C++ equivalent.
	MCSampler* SamplerConstruct(const mxArray* structptr, Hilbert* hlb)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Nsamp";
		matvar = FindField(structptr, fieldname);
		int Nsamp = mex2int(matvar);
		fieldname = "Nequil";
		matvar = FindField(structptr, fieldname);
		int Nequil = mex2int(matvar);
		fieldname = "Nblock";
		matvar = FindField(structptr, fieldname);
		int Nblock = mex2int(matvar);
		fieldname = "Hamiltonian";
		matvar = FindField(structptr, fieldname);
		Hamiltonian* hmt = HamiltonianConstruct(matvar);
		MCSampler* mcs = new MCSampler(hlb, hmt, Nsamp, Nequil, Nblock);
		return mcs;
	}
	// SamplerExtract will search the file for the given struct name and attempt to construct a MCSampler object.
	// - Requires a Hilbert instance to perform.
	MCSampler* SamplerExtract(const char* filename, const char* structname, Hilbert* hlb)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		MCSampler* mcs = SamplerConstruct(propstruct, hlb);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (mcs == NULL)
		{
			cerr << "Unable to construct appropriate Ansatz object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New MCSampler constructed successfully." << endl;
		return mcs;
	}

	// SRConstruct takes a SR Matlab struct and outputs the C++ equivalent.
	StochReconfig* SRConstruct(const mxArray* structptr)
	{
		const mxArray* matvar;
		const char* fieldname;
		fieldname = "Npass";
		matvar = FindField(structptr, fieldname);
		int Npass = mex2int(matvar);
		fieldname = "Ncore";
		matvar = FindField(structptr, fieldname);
		int Ncore = mex2int(matvar);
		fieldname = "ExtraSamp";
		matvar = FindField(structptr, fieldname);
		int ExtraSamp = mex2int(matvar);
		fieldname = "LearnRate";
		matvar = FindField(structptr, fieldname);
		double LearnRate = mex2dbl(matvar);
		fieldname = "Regularisation";
		matvar = FindField(structptr, fieldname);
		vector<double> Regularisation = mex2std_dvec(matvar);
		fieldname = "Rollback";
		matvar = FindField(structptr, fieldname);
		vector<int> Rollback = mex2std_ivec(matvar);
		fieldname = "EnergyTolerance";
		matvar = FindField(structptr, fieldname);
		vector<double> EnergyTolerance = mex2std_dvec(matvar);
		fieldname = "EnergySensitivity";
		matvar = FindField(structptr, fieldname);
		vector<double> EnergySensitivity = mex2std_dvec(matvar);
		fieldname = "ParamSensitivity";
		matvar = FindField(structptr, fieldname);
		vector<double> ParamSensitivity = mex2std_dvec(matvar);
		StochReconfig* sr_ = new StochReconfig(Npass, Ncore, ExtraSamp, LearnRate, 
			Regularisation, Rollback, EnergyTolerance, EnergySensitivity, ParamSensitivity);
		return sr_;
	}
	// SRExtract will search the file for the given struct name and attempt to construct a StochReconfig object.
	StochReconfig* SRExtract(const char* filename, const char* structname)
	{
		MATFile* matfile = matOpen(filename, "r");
		if (matfile == NULL)
		{
			cerr << "Error reading file " << filename << ". Check name is correct." << endl;
			std::abort();
		}
		const mxArray* propstruct;
		propstruct = matGetVariable(matfile, structname);
		if (!mxIsClass(propstruct, "struct"))
		{
			cerr << "Target variable is not a struct." << endl;
			std::abort();
		}
		StochReconfig* sr_ = SRConstruct(propstruct);
		if (matClose(matfile) != 0)
		{
			cerr << "Error closing file." << endl;
			std::abort();
		}
		if (sr_ == NULL)
		{
			cerr << "Unable to construct appropriate Ansatz object - check struct fields match." << endl;
			std::abort();
		}
		cout << "New StochReconfig constructed successfully." << endl;
		return sr_;
	}

	// VectorMatSave will take a vector (Eig or std) and save it to the given filename with the given variable name.
	// - Eigen vector input version:
	void VectorMatSave(VectorXd EigVec, const char* filename, const char* varname)
	{
		MATFile* savemat = matOpen(filename, "u");
		if (savemat == NULL)
		{
			cout << "No file with the given file name exists - writing new file..." << endl;
			savemat = matOpen(filename, "w");
			if (savemat == NULL)
			{
				cerr << "Error creating file " << filename << " - check permissions." << endl;
				std::abort();
			}
		}
		size_t lv = EigVec.size();
		mxArray* savemx = mxCreateDoubleMatrix(lv, 1, mxREAL);
		if (savemx == NULL)
		{
			cerr << "Error creating mxArray " << varname << " - check memory is sufficient." << endl;
			std::abort();
		}
		memcpy(mxGetPr(savemx), EigVec.data(), lv * sizeof(double));
		matError status = matPutVariable(savemat, varname, savemx);
		if (status != 0)
		{
			cerr << "Error writing variable " << varname << " into file " << filename << "." << endl;
			std::abort();
		}
		else
		{
			cout << "Variable " << varname << " written into file " << filename << "." << endl;
		}
		if (matClose(savemat) != 0)
		{
			cerr << "Error closing .mat file." << endl;
			std::abort();
		}
		return;
	}
	// - std vector input version:
	void VectorMatSave(vector<double> StdVec, const char* filename, const char* varname)
	{
		MATFile* savemat = matOpen(filename, "u");
		if (savemat == NULL)
		{
			cout << "No file with the given file name exists - writing new file..." << endl;
			savemat = matOpen(filename, "w");
			if (savemat == NULL)
			{
				cerr << "Error creating file " << filename << " - check permissions." << endl;
				std::abort();
			}
		}
		size_t lv = StdVec.size();
		mxArray* savemx = mxCreateDoubleMatrix(lv, 1, mxREAL);
		if (savemx == NULL)
		{
			cerr << "Error creating mxArray " << varname << " - check memory is sufficient." << endl;
			std::abort();
		}
		memcpy(mxGetPr(savemx), StdVec.data(), lv * sizeof(double));
		matError status = matPutVariable(savemat, varname, savemx);
		if (status != 0)
		{
			cerr << "Error writing variable " << varname << " into file " << filename << "." << endl;
			std::abort();
		}
		else
		{
			cout << "Variable " << varname << " written into file " << filename << "." << endl;
		}
		if (matClose(savemat) != 0)
		{
			cerr << "Error closing .mat file." << endl;
			std::abort();
		}
		return;
	}

	// MatrixMatSave will take a matrix (Eig or std) and save it to the given filename with the given variable name.
	// - Eigen matrix input version:
	void MatrixMatSave(MatrixXd EigMat, const char* filename, const char* varname)
	{
		MATFile* savemat = matOpen(filename, "u");
		if (savemat == NULL)
		{
			cout << "No file with the given file name exists - writing new file..." << endl;
			savemat = matOpen(filename, "w");
			if (savemat == NULL)
			{
				cerr << "Error creating file " << filename << " - check permissions." << endl;
				std::abort();
			}
		}
		mxArray* savemx = mxCreateDoubleMatrix(EigMat.rows(), EigMat.cols(), mxREAL);
		if (savemx == NULL)
		{
			cerr << "Error creating mxArray " << varname << " - check memory is sufficient." << endl;
			std::abort();
		}
		memcpy(mxGetPr(savemx), EigMat.data(), EigMat.size() * sizeof(double));
		matError status = matPutVariable(savemat, varname, savemx);
		if (status != 0)
		{
			cerr << "Error writing variable " << varname << " into file " << filename << "." << endl;
			std::abort();
		}
		else
		{
			cout << "Variable " << varname << " written into file " << filename << "." << endl;
		}
		if (matClose(savemat) != 0)
		{
			cerr << "Error closing .mat file." << endl;
			std::abort();
		}
		return;
	}
	// - std vector input version:
	void MatrixMatSave(vector<vector<double>> StdMat, const char* filename, const char* varname)
	{
		MATFile* savemat = matOpen(filename, "u");
		if (savemat == NULL)
		{
			cout << "No file with the given file name exists - writing new file..." << endl;
			savemat = matOpen(filename, "w");
			if (savemat == NULL)
			{
				cerr << "Error creating file " << filename << " - check permissions." << endl;
				std::abort();
			}
		}
		// Assume column major i.e. each nested vector is a column.
		mxArray* savemx = mxCreateDoubleMatrix(StdMat[0].size(), StdMat.size(), mxREAL);
		if (savemx == NULL)
		{
			cerr << "Error creating mxArray " << varname << " - check memory is sufficient." << endl;
			std::abort();
		}
		// For ease of use, assign to a temporary vector then extract data.
		size_t numel = StdMat.size() * StdMat[0].size();
		vector<double> StdVec(numel);
		int m = 0;
		for (int j = 0; j < StdMat.size(); j++)
		{
			for (int i = 0; i < StdMat[j].size(); i++)
			{
				StdVec[m] = StdMat[j][i];
				m++;
			}
		}
		memcpy(mxGetPr(savemx), StdVec.data(), numel * sizeof(double));
		matError status = matPutVariable(savemat, varname, savemx);
		if (status != 0)
		{
			cerr << "Error writing variable " << varname << " into file " << filename << "." << endl;
			std::abort();
		}
		else
		{
			cout << "Variable " << varname << " written into file " << filename << "." << endl;
		}
		if (matClose(savemat) != 0)
		{
			cerr << "Error closing .mat file." << endl;
			std::abort();
		}
		return;
	}

	// ScalarMatSave will take a double input and save/append it to the given filename with the given variable name.
	void ScalarMatSave(double StdDbl, const char* filename, const char* varname)
	{
		MATFile* savemat = matOpen(filename, "u");
		if (savemat == NULL)
		{
			cout << "No file with the given file name exists - writing new file..." << endl;
			savemat = matOpen(filename, "w");
			if (savemat == NULL)
			{
				cerr << "Error creating file " << filename << " - check permissions." << endl;
				std::abort();
			}
		}
		// Assume column major i.e. each nested vector is a column.
		mxArray* savemx = mxCreateDoubleMatrix(1, 1, mxREAL);
		if (savemx == NULL)
		{
			cerr << "Error creating mxArray " << varname << " - check memory is sufficient." << endl;
			std::abort();
		}
		// For ease of use, assign to a temporary vector then extract data.
		vector<double> StdVec(1);		
		StdVec[0] = StdDbl;
		memcpy(mxGetPr(savemx), StdVec.data(), sizeof(double));
		matError status = matPutVariable(savemat, varname, savemx);
		if (status != 0)
		{
			cerr << "Error writing variable " << varname << " into file " << filename << "." << endl;
			std::abort();
		}
		else
		{
			cout << "Variable " << varname << " written into file " << filename << "." << endl;
		}
		if (matClose(savemat) != 0)
		{
			cerr << "Error closing .mat file." << endl;
			std::abort();
		}
		return;
	}
}

#endif