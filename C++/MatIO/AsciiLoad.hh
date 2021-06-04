// ------------------------------------------------------------------------------------------------------------
// AsciiLoad.hh - contains the processes necessary to convert from loaded ASCII files containing parameters
//				  to the C++ Eigen vectors and back. Trial for interfacing with the Matlab code.
// ------------------------------------------------------------------------------------------------------------

#ifndef ASCIILOAD
#define ASCIILOAD

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

namespace nqsvmc
{
	VectorXd ParamLoad_ASCII(const char* filename)
	{
		ifstream file(filename);
		double val_;
		vector<double> val_vec;
		if (file.is_open())
		{
			while (file >> val_)
			{
				val_vec.push_back(val_);
			}
			file.close();
		}
		Map<VectorXd> val_eig(&val_vec[0], val_vec.size());
		return val_eig;
	}

	void ParamSave_ASCII(VectorXd Params, const char* filename)
	{
		ofstream file(filename, std::ofstream::out);
		if (file.is_open())
		{
			for (int p = 0; p < Params.size(); p++)
			{
				file << Params(p) << "\n";
			}
		}
		file.close();
		std::cout << "Parameter vector saved as " << filename << ".\n";
		return;
	}

	void ParamSave_ASCII(MatrixXd Params, const char* filename)
	{
		ofstream file(filename, std::ofstream::out);
		if (file.is_open())
		{
			for (int p = 0; p < Params.size(); p++)
			{
				file << Params(p) << "\n";
			}
		}
		file.close();
		std::cout << "Parameter vector saved as " << filename << ".\n";
		return;
	}
}

#endif