// ------------------------------------------------------------------------------------------------------------
// Graph.cpp - contains the definitions of the Graph object class, which involves a wrapper implementer 
//			   subclass as an interface for the rest of the code.
// ------------------------------------------------------------------------------------------------------------

#ifndef GRAPH_CPP
#define GRAPH_CPP

#include <Eigen/Dense>
#include "Graph.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// Constructors for the Graph wrapper class:
	// - Generate a HypCub subclass instance.
	Graph::Graph(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound, char ID)
	{
		if (ID == 'c') // 'c' denotes a hypercube
		{
			g_ = new HypCub(dim, lvecs, bound);
		}
		// To add - honeycomb (h), triangular (t).
	}
	// - Custom construction with premade lookup lists.
	Graph::Graph(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound, vector<vector<int>> bonds,
		vector<int> slinds, bool bondmapflag)
	{
		// If provided the full set of information, generates a CustomGraph.
		g_ = new CustomGraph(dim, lvecs, bound, bonds, slinds, bondmapflag);
	}
	// Observer functions for extracting information from a general Graph object
	// - Nsite will return the number of sites in the lattice.
	int Graph::Nsite() const
	{
		return g_->Nsite();
	}
	// - GraphDim will return the dimensions defined in the Graph.
	vector<int> Graph::GraphDim() const
	{
		return g_->GraphDim();
	}
	// - Ntranslate will return the number of lookup lists defined in Graph for different translations.
	int Graph::Ntranslate() const
	{
		return g_->Ntranslate();
	}
	// - Nvecs will return the number of primitive lattice vectors defined in the Graph.
	int Graph::Nvecs() const
	{
		return g_->Nvecs();
	}
	// - Coord will return the coordination of the lattice.
	int Graph::Coord() const
	{
		return g_->Coord();
	}
	// - LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
	vector<int> Graph::LatVec(int index) const
	{
		return g_->LatVec(index);
	}
	// - VecInd returns the lattice vector coefficients associated with the input translation index.
	vector<int> Graph::VecInd(int index) const
	{
		return g_->VecInd(index);
	}
	// - BondSearch returns the lookup list associated with the input translation index.
	vector<int> Graph::BondSearch(int index) const
	{
		return g_->BondSearch(index);
	}
	// - BondRead returns the primary neighbour lookup list associated with the Graph.
	vector<vector<int>> Graph::BondRead() const
	{
		return g_->BondRead();
	}
	// - SLRead returns the stored sublattice indices associated with the sites in Graph.
	vector<int> Graph::SLRead() const
	{
		return g_->SLRead();
	}
	// - SLLoad alters the existing sublattice indices associated with the sites in Graph.
	void Graph::SLLoad(vector<int> indices)
	{
		return g_->SLLoad(indices);
	}
	// - Graph2Matrix outputs an Eigen matrix of connections using primary bonds.
	MatrixXd Graph::Graph2Matrix() const
	{
		int N = g_->Nsite();
		int Nb = g_->Nvecs();
		MatrixXd ConMat = MatrixXd::Zero(N, N);
		vector<vector<int>> Bonds = g_->BondRead();
		for (int b = 0; b < Nb; b++)
		{
			for (int n = 0; n < N; n++)
			{
				int Site = Bonds[b][n];
				if (Site >= 0)
				{
					ConMat(n, Site) = 1;
					ConMat(Site, n) = 1;
				}
			}
		}
		return ConMat;
	}

	// <<<<<<------------ GeneralGraph functions ------------>>>>>>

	// BondGen outputs a nested vector list of site links, format is link = [vector index][site]
	vector<vector<int>> GeneralGraph::BondGen(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound)
	{
		vector<vector<int>> Bonds;
		int N = 1;
		vector<int> dimprod; // Cumulative product of dimensions.
		dimprod.resize(dim.size());
		int x, y, dx;
		for (int d = 0; d < dim.size(); d++)
		{
			dimprod[d] = N;
			N *= dim[d];
		}
		Bonds.resize(lvecs.size()); // One lookup list per provided vector.
		for (int v = 0; v < lvecs.size(); v++)
		{
			Bonds[v].resize(N); // Ensure all lookup lists are correct size.
			for (int n = 0; n < N; n++)
			{
				Bonds[v][n] = n;
			}
		}
		// Begin populating lists.
		for (int d = 0; d < dim.size(); d++)
		{
			for (int v = 0; v < lvecs.size(); v++)
			{
				for (int n = 0; n < N; n++)
				{
					// For each entry, establish co-ordinate in current dimension.
					if (Bonds[v][n] >= 0)
					{
						x = (int)floor(Bonds[v][n] / dimprod[d]) % dim[d];
						y = x + lvecs[v][d];
						if (bound[d] == 0)
						{
							if ((y >= dim[d]) || (y < 0))
							{
								Bonds[v][n] = -1; // Treat negative entries as empty / no connection.
							}
							else
							{
								dx = (y - x) * dimprod[d];
								Bonds[v][n] += dx;
							}
						}
						else
						{
							y = y % dim[d];
							dx = (y - x) * dimprod[d];
							Bonds[v][n] += dx;
						}
					}
				}
			}
		}
		return Bonds;
	}

	// BondMapGen populates the BondMap field of a general Graph, for a full connectivity list.
	void GeneralGraph::BondMapGen(vector<int> dim, vector<vector<int>> lvecs, vector<vector<int>> bonds,
		vector<vector<int>>* vecinds, vector<vector<int>>* bondmap)
	{ // Given some principal lattice vectors and associated lookup tables, will create a full lookup list.
	  // Assumes each of the entries in bonds corresponds to a lattice vector.
		int nvec = (int)lvecs.size();
		int N = 1;
		for (int d = 0; d < dim.size(); d++)
		{
			N *= dim[d];
		}
		int Ntr = 1;
		vector<int> vtr(nvec);
		vector<int> ptr(nvec); // Used for indexing the various combinations of vectors.
		vector<int> tvec(dim.size());
		for (int v = 0; v < nvec; v++)
		{
			ptr[v] = Ntr;
			tvec = lvecs[v];
			int s = 1; // Work out number of translations required to map onto same site with PBC.
			int rem = 0;
			for (int d = 0; d < dim.size(); d++)
			{
				rem += lvecs[v][d];
			}
			while (rem > 0)
			{
				s++;
				int trem = 0;
				for (int d = 0; d < dim.size(); d++)
				{
					tvec[d] += lvecs[v][d];
					trem += tvec[d] % dim[d];
				}
				rem = trem;
			}
			vtr[v] = s;
			Ntr *= vtr[v];
		}
		vector<bool> vlist(Ntr); // Vector listing valid vector combinations to avoid double counting.
		vector<vector<int>> KOList(N, vector<int>(N)); // Knockout list to avoid double counting mappings.
		vector<vector<int>> vecindsT(Ntr, vector<int>(nvec)); // Temporary vector index.
		vector<vector<int>> bondmapT(Ntr, vector<int>(N)); // Temporary bond map.
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < N; m++)
			{
				KOList[n][m] = m;
			}
		}
		int dest; // Destination site and map index.
		int NtrA = Ntr;
		bool delflag = false;
		// Start combining lattice vectors and creating the lookup lists.
		for (int n = 0; n < Ntr; n++)
		{
			for (int v = 0; v < lvecs.size(); v++)
			{
				vecindsT[n][v] = (int)floor(n / ptr[v]);
			}
			for (int m = 0; m < N; m++)
			{
				dest = m;
				for (int x = 0; x < lvecs.size(); x++)
				{
					for (int dx = 0; dx < vecindsT[n][x]; dx++)
					{
						dest = bonds[x][dest];
						if (dest < 0)
						{
							break;
						}
					}
				}
				bondmapT[n][m] = dest;
				if (dest >= 0)
				{
					if (KOList[m][dest] < 0)
					{
						delflag = true;
					}
					else
					{
						KOList[m][dest] = -1;
						delflag = false;
					}
				}
			}
			vlist[n] = delflag;
			NtrA -= delflag;
		}
		// Trim redundant / excess translation combinations and log lattice vector combinations.
		vecinds->resize(NtrA);
		bondmap->resize(NtrA);
		Ntr = NtrA;
		int c = 0;
		for (int m = 0; m < Ntr; m++)
		{
			if (!vlist[m])
			{
				vecinds[0][c] = vecindsT[m];
				bondmap[0][c] = bondmapT[m];
				c++;
			}
		}
		return;
	}

	

	// Definitions for Graph subclass functions here:

	// <<<<<<------------ HypCub functions ------------>>>>>>

	// - HypCub constructor with all inputs (dimensions, vectors and boundary condition):
	HypCub::HypCub(vector<int> dimensions, vector<vector<int>> vectors, vector<int> boundcon)
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
	// - HypCub::Nsite will return the number of sites in the lattice.
	int HypCub::Nsite() const
	{
		return N;
	}
	// - HypCub::GraphDim will return the dimensions defined in the Graph.
	vector<int> HypCub::GraphDim() const
	{
		return Dim;
	}
	// - HypCub::Ntranslate will return the number of lookup lists defined in Graph for different translations.
	int HypCub::Ntranslate() const
	{
		return Ntr;
	}
	// - HypCub::Nvecs will return the number of primitive lattice vectors defined in the Graph.
	int HypCub::Nvecs() const
	{
		return (int)LVecs.size();
	}
	// - HypCub::Coord will return the coordination number of the lattice.
	int HypCub::Coord() const
	{
		return z;
	}
	// - HypCub::LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
	vector<int> HypCub::LatVec(int index) const
	{
		if ((index >= LVecs.size()) || (index < 0))
		{
			cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
			std::abort();
		}
		return LVecs[index];
	}
	// - HypCub::VecInd returns the lattice vector coefficients associated with the input translation index.
	vector<int> HypCub::VecInd(int index) const
	{
		if ((index >= LVecs.size()) || (index < 0))
		{
			cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
			std::abort();
		}
		return VInds[index];
	}
	// - HypCub::BondMap returns the lookup list associated with the input translation index.
	vector<int> HypCub::BondSearch(int index) const
	{
		if ((index > Ntr) || (index < 0))
		{
			cerr << "Invalid translation index. Use Ntr to determine the number of distinct translations." << endl;
			std::abort();
		}
		return BondMap[index];
	}
	// - HypCub::BondRead returns the primary neighbour lookup list associated with the Graph.
	vector<vector<int>> HypCub::BondRead() const
	{
		return Bonds;
	}
	// - HypCub::SLRead returns the stored sublattice indices associated with the sites in Graph.
	vector<int> HypCub::SLRead() const
	{
		return SLInds;
	}
	// - HypCub::SLLoad alters the existing sublattice indices associated with the sites in Graph.
	void HypCub::SLLoad(vector<int> indices)
	{
		if (indices.size() != N)
		{
			cerr << "Input vector does not have the correct number of entries." << endl;
			std::abort();
		}
		swap(indices, SLInds);
	}

	// <<<<<<------------ CustomGraph functions ------------>>>>>>
	// - Constructor with full inputs (dimensions, vectors and boundary condition):
	CustomGraph::CustomGraph(vector<int> dimensions, vector<vector<int>> vectors, vector<int> boundcon, vector<vector<int>> bonds,
		vector<int> sublatinds, bool bondmapflag)
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
		Bonds = bonds;
		if (bondmapflag)
		{
			BondMapGen(dimensions, vectors, Bonds, &VInds, &BondMap);
		}
		Ntr = (int)BondMap.size();
		SLInds = sublatinds;
	}
	// - CustomGraph::Nsite will return the number of sites in the lattice.
	int CustomGraph::Nsite() const
	{
		return N;
	}
	// - CustomGraph::GraphDim will return the dimensions defined in the Graph.
	vector<int> CustomGraph::GraphDim() const
	{
		return Dim;
	}
	// - CustomGraph::Ntranslate will return the number of lookup lists defined in Graph for different translations.
	int CustomGraph::Ntranslate() const
	{
		return Ntr;
	}
	// - CustomGraph::Nvecs will return the number of primitive lattice vectors defined in the Graph.
	int CustomGraph::Nvecs() const
	{
		return (int)LVecs.size();
	}
	// - CustomGraph::Coord will return the coordination number of the lattice.
	int CustomGraph::Coord() const
	{
		return z;
	}
	// - CustomGraph::LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
	vector<int> CustomGraph::LatVec(int index) const
	{
		if ((index >= LVecs.size()) || (index < 0))
		{
			cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
			std::abort();
		}
		return LVecs[index];
	}
	// - CustomGraph::VecInd returns the lattice vector coefficients associated with the input translation index.
	vector<int> CustomGraph::VecInd(int index) const
	{
		if ((index >= LVecs.size()) || (index < 0))
		{
			cerr << "Invalid vector index. Use Nvecs to determine the number of principal lattice vectors." << endl;
			std::abort();
		}
		return VInds[index];
	}
	// - CustomGraph::BondMap returns the lookup list associated with the input translation index.
	vector<int> CustomGraph::BondSearch(int index) const
	{
		if ((index > Ntr) || (index < 0))
		{
			cerr << "Invalid translation index. Use Ntr to determine the number of distinct translations." << endl;
			std::abort();
		}
		return BondMap[index];
	}
	// - CustomGraph::BondRead returns the primary neighbour lookup list associated with the Graph.
	vector<vector<int>> CustomGraph::BondRead() const
	{
		return Bonds;
	}
	// - CustomGraph::SLRead returns the stored sublattice indices associated with the sites in Graph.
	vector<int> CustomGraph::SLRead() const
	{
		return SLInds;
	}
	// - CustomGraph::SLLoad alters the existing sublattice indices associated with the sites in Graph.
	void CustomGraph::SLLoad(vector<int> indices)
	{
		if (indices.size() != N)
		{
			cerr << "Input vector does not have the correct number of entries." << endl;
			std::abort();
		}
		swap(indices, SLInds);
	}
}

#endif