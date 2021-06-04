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
	// Defining some GeneralGraph helper functions:
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
}

#endif