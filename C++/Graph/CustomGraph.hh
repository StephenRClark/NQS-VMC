// ------------------------------------------------------------------------------------------------------------
// CustomGraph.hh - contains the definitions of the CustomGraph subclass, which is generated from pre-made 
//					look-up lists and bundled together. These are commonly made when loading in data.
// ------------------------------------------------------------------------------------------------------------

#ifndef CUSTOMGRAPH_HH
#define CUSTOMGRAPH_HH

#include <iostream>
#include <vector>
#include "Graph.hh"

using namespace std;

namespace nqsvmc
{
	class CustomGraph : public GeneralGraph // Basic hypercube Graph.
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
		// - Constructor with full inputs (dimensions, vectors and boundary condition):
		CustomGraph(vector<int> dimensions, vector<vector<int>> vectors, vector<int> boundcon, vector<vector<int>> bonds,
			vector<int> sublatinds, bool bondmapflag);
		// Observer functions for extracting information from a general Graph object
		// - Nsite will return the number of sites in the lattice.
		int Nsite() const;
		// - GraphDim will return the dimensions defined in the Graph.
		vector<int> GraphDim() const;
		// - Ntranslate will return the number of lookup lists defined in Graph for different translations.
		int Ntranslate() const;
		// - Nvecs will return the number of primitive lattice vectors defined in the Graph.
		int Nvecs() const;
		// - Coord will return the coordination number of the lattice.
		int Coord() const;
		// - LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
		vector<int> LatVec(int index) const;
		// - VecInd returns the lattice vector coefficients associated with the input translation index.
		vector<int> VecInd(int index) const;
		// - BondMap returns the lookup list associated with the input translation index.
		vector<int> BondSearch(int index) const;
		// - BondRead returns the primary neighbour lookup list associated with the Graph.
		vector<vector<int>> BondRead() const;
		// - SLRead returns the stored sublattice indices associated with the sites in Graph.
		vector<int> SLRead() const;
		// - SLLoad alters the existing sublattice indices associated with the sites in Graph.
		void SLLoad(vector<int> indices);
	};
}

#endif