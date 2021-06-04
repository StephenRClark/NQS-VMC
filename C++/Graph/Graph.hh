// ------------------------------------------------------------------------------------------------------------
// Graph.hh - contains the definitions of the Graph object class, which represents the dimensions and 
//			  connectivity of the lattice. Basically acts as a collection of lookup lists.
// ------------------------------------------------------------------------------------------------------------

#ifndef GRAPH_HH
#define GRAPH_HH

#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class GeneralGraph // Abstract layout for the framework of the object class.
	{
	public:
		// Observer functions for extracting information from a general Graph object
		// - Nsite will return the number of sites in the lattice.
		virtual int Nsite() const = 0;
		// - GraphDim will return the dimensions defined in the Graph.
		virtual vector<int> GraphDim() const = 0;
		// - Ntranslate will return the number of lookup lists defined in Graph for different translations.
		virtual int Ntranslate() const = 0;
		// - Nvecs will return the number of primitive lattice vectors defined in the Graph.
		virtual int Nvecs() const = 0;
		// - Coord will return the coordination number of the lattice.
		virtual int Coord() const = 0;
		// - LatVec returns one of the lattice vectors defined in the Graph, depending on the input index value.
		virtual vector<int> LatVec(int index) const = 0;
		// - VecInd returns the lattice vector coefficients associated with the input translation index.
		virtual vector<int> VecInd(int index) const = 0;
		// - BondSearch returns the lookup list associated with the input translation index.
		virtual vector<int> BondSearch(int index) const = 0;
		// - BondRead returns the primary neighbour lookup list associated with the Graph.
		virtual vector<vector<int>> BondRead() const = 0;
		// - SLRead returns the stored sublattice indices associated with the sites in Graph.
		virtual vector<int> SLRead() const = 0;
		// - SLLoad alters the existing sublattice indices associated with the sites in Graph.
		virtual void SLLoad(vector<int> indices) = 0;

		// Forward declaration of some Graph helper functions.
		vector<vector<int>> BondGen(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound);
		void BondMapGen(vector<int> dim, vector<vector<int>> lvecs, vector<vector<int>> bonds,
			vector<vector<int>>* vecinds, vector<vector<int>>* bondmap);
	};

	// Forward declaration of implementer subclass:
	class Graph : public GeneralGraph
	{
	private:
		GeneralGraph* g_;
	public:
		// Constructors for various Graph subclasses.
		// - Preloaded subclass constructor - supply identifier character ID to select.
		Graph(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound, char ID);
		// - Custom construction with premade lookup lists.
		Graph(vector<int> dim, vector<vector<int>> lvecs, vector<int> bound, vector<vector<int>> bonds,
			vector<int> slinds, bool bondmapflag);
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
		// - BondSearch returns the lookup list associated with the input translation index.
		vector<int> BondSearch(int index) const;
		// - BondRead returns the primary neighbour lookup list associated with the Graph.
		vector<vector<int>> BondRead() const;
		// - SLRead returns the stored sublattice indices associated with the sites in Graph.
		vector<int> SLRead() const;
		// - SLLoad alters the existing sublattice indices associated with the sites in Graph.
		void SLLoad(vector<int> indices);
		// - Graph2Matrix outputs an Eigen matrix of connections using primary bonds.
		MatrixXd Graph2Matrix() const;
	};
}

#include "HypCub.hh"
#include "CustomGraph.hh"

#endif
