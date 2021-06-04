// ------------------------------------------------------------------------------------------------------------
// HypCub.hh - contains the definitions of commonly used Graph subclasses that have well known structure.
//			   Expect to eventually add more to this list.
// ------------------------------------------------------------------------------------------------------------

#ifndef HYPCUB_HH
#define HYPCUB_HH

#include <iostream>
#include <vector>
#include "Graph.hh"

using namespace std;

namespace nqsvmc
{
	class HypCub : public GeneralGraph // Basic hypercube Graph.
	{
	private:
		vector<int> Dim; // Dimensions of the lattice.
		int N; // Total number of sites.
		vector<vector<int>> LVecs; // Primary lattice vectors.
		int z; // Co-ordination number.
		vector<vector<int>> Bonds; // Primary neighbour list - z/2 entries (no double counting of bonds).
		vector<int> Bound; // Boundary conditions along each dimension. 0 = closed, 1 = periodic.
		int Ntr; // Number of unique translations possible with the given primary lattice vectors.
		vector<vector<int>> VInds; // List of lattice vector coefficients associated with a translation index - Ntr entries.
		vector<vector<int>> BondMap; // Nested vector array of lookup lists associated with a translation index - Ntr entries.
		vector<int> SLInds; // List of sublattice indices - N entries.
	public:
		// Constructor functions:
		// - Constructor with proper inputs (dimensions, vectors and boundary condition):
		HypCub(vector<int> dimensions, vector<vector<int>> vectors, vector<int> boundcon)
		{
			Dim = dimensions;
			N = 1;
			for (int d = 0; d < dimensions.size(); d++)
			{
				N *= dimensions[d];
			}
			z = 2 * (int)dimensions.size();
			LVecs = vectors;
			if (dimensions.size() != boundcon.size())
			{
				cerr << "Number of boundary conditions does not match number of dimensions." << endl;
				std::abort();
			}
			Bound = boundcon;
			Bonds = BondGen(dimensions, vectors, boundcon);
			BondMapGen(dimensions, vectors, Bonds, &VInds, &BondMap);
			Ntr = (int)BondMap.size();
			SLInds.resize(N);
			for (int n = 0; n < N; n++)
			{
				SLInds[n] = 1;
			}
		}
		// Observer functions for extracting information from a general Graph object
		// - Nsite will return the number of sites in the lattice.
		int Nsite() const
		{
			return N;
		}
		// - GraphDim will return the dimensions defined in the Graph.
		vector<int> GraphDim() const
		{
			return Dim;
		}
		// - Ntranslate will return the number of lookup lists defined in Graph for different translations.
		int Ntranslate() const
		{
			return Ntr;
		}
		// - Nvecs will return the number of primitive lattice vectors defined in the Graph.
		int Nvecs() const
		{
			return (int)LVecs.size();
		}
		// - Coord will return the coordination number of the lattice.
		int Coord() const
		{
			return z;
		}
		// - LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
		vector<int> LatVec(int index) const
		{
			if ((index >= LVecs.size()) || (index < 0))
			{
				cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
				std::abort();
			}
			return LVecs[index];
		}
		// - VecInd returns the lattice vector coefficients associated with the input translation index.
		vector<int> VecInd(int index) const
		{
			if ((index >= LVecs.size()) || (index < 0))
			{
				cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
				std::abort();
			}
			return VInds[index];
		}
		// - BondMap returns the lookup list associated with the input translation index.
		vector<int> BondSearch(int index) const
		{
			if ((index > Ntr) || (index < 0))
			{
				cerr << "Invalid translation index. Use Ntr to determine the number of distinct translations." << endl;
				std::abort();
			}
			return BondMap[index];
		}
		// - Bonds returns the primary neighbour lookup list associated with the Graph.
		vector<vector<int>> BondRead() const
		{
			return Bonds;
		}
		// - SLRead returns the stored sublattice indices associated with the sites in Graph.
		vector<int> SLRead() const
		{
			return SLInds;
		}
		// - SLLoad alters the existing sublattice indices associated with the sites in Graph.
		void SLLoad(vector<int> indices)
		{
			if (indices.size() != N)
			{
				cerr << "Input vector does not have the correct number of entries." << endl;
				std::abort();
			}
			swap(indices, SLInds);
		}
	};
}

#endif