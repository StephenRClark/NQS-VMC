// ------------------------------------------------------------------------------------------------------------
// nqsvmc.cpp - the main source file, pulling in all the necessary definitions from the headers and executing 
//				the actual VMC process.
// ------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <chrono>
#include <string>
#include <sstream>
#include <random>
#include <Eigen/Dense>
#include "nqsvmc.hh"

using namespace std;
using namespace std::chrono;
using namespace Eigen;
using namespace nqsvmc;

int main(int argc, char* argv[])
{	
	// Command line arguments: 1) name of a text file containing setup filenames,
	//						   2) name of a textfile containing the wavefunction filenames.
	if (argc < 3)
	{
		cerr << "Please specify both the setup and wavefunction filename lists." << endl;
		abort();
	}
	if (argc > 3)
	{
		cerr << "Too many files specified - if setting up multiple runs, create text files listing the setup and wavefunction filenames." << endl;
		abort();
	}
	// Load in setup text file and look at number of files that need loading.
	ifstream setuplist(argv[1]);
	vector<string> setupfilenames;
	string filename;
	if (setuplist.is_open())
	{
		while (std::getline(setuplist, filename))
		{
			setupfilenames.push_back(filename);
		}
		setuplist.close();
	}
	size_t Nsetup = setupfilenames.size();
	cout << Nsetup << " setup filenames found." << endl;
	ifstream wfunclist(argv[2]);
	vector<string> wfuncfilenames;
	if (wfunclist.is_open())
	{
		while (std::getline(wfunclist, filename))
		{
			wfuncfilenames.push_back(filename);
		}
		wfunclist.close();
	}
	size_t Nwfunc = wfuncfilenames.size();
	cout << Nwfunc << " wavefunction filenames found." << endl;
	size_t Nrun = Nwfunc;
	if (Nsetup != Nwfunc)
	{
		cout << "Setup file and wavefunction file number mismatch." << endl;
		Nrun = min(Nwfunc, Nsetup);
		cout << "Only " << Nrun << " optimisations can be performed." << endl;
	}
	cout << "Current implementation assumes a three-stage optimisation process - check setup files in case of errors." << endl;
	Hilbert* hlb;
	MCSampler* Samp1;
	MCSampler* Samp2;
	MCSampler* Samp3;
	MCSampler* EvalSamp;
	StochReconfig* SR1;
	StochReconfig* SR2;
	StochReconfig* SR3;
	Ansatz* Ans;
	VectorXd EnIter;
	vector<Operator*> EvalOps;
	vector<const char*> OpNames;
	vector<const char*> ValNames;
	OpNames.push_back("EnSqProp"); ValNames.push_back("EnSq");
	OpNames.push_back("VarNProp"); ValNames.push_back("VarN");
	OpNames.push_back("NiNjProp"); ValNames.push_back("NiNj");
	OpNames.push_back("DbHlProp"); ValNames.push_back("DbHl");
	OpNames.push_back("OcFrProp"); ValNames.push_back("OcFr");
	OpNames.push_back("BiBjProp"); ValNames.push_back("BiBj");
	vector<MatrixXd> EvalVals;
	double EneGS;
	double RunTime;
	double EvalTime;
	high_resolution_clock::time_point tstart;
	high_resolution_clock::time_point tfinish;
	duration<double> trun;
	for (int r = 0; r < Nrun; r++)
	{
		std::cout << endl << "Performing optimisation of " << wfuncfilenames[r] << " with " << setupfilenames[r] << "." << endl << endl;
		const char* varname;
		std::cout << "Loading Hilbert object..." << endl;
		varname = "HilbertProp";
		hlb = HilbertExtract(setupfilenames[r].c_str(), varname);
		std::cout << "Loading Sampler objects... " << endl;
		varname = "Samp1Prop";
		std::cout << "1... ";
		Samp1 = SamplerExtract(setupfilenames[r].c_str(), varname, hlb);
		varname = "Samp2Prop";
		std::cout << "2... ";
		Samp2 = SamplerExtract(setupfilenames[r].c_str(), varname, hlb);
		std::cout << "3... ";
		varname = "Samp3Prop";
		Samp3 = SamplerExtract(setupfilenames[r].c_str(), varname, hlb);
		std::cout << endl << "Loading SR objects... " << endl;
		varname = "SR1Prop";
		std::cout << "1... ";
		SR1 = SRExtract(setupfilenames[r].c_str(), varname);
		varname = "SR2Prop";
		std::cout << "2... ";
		SR2 = SRExtract(setupfilenames[r].c_str(), varname);
		varname = "SR3Prop";
		std::cout << "3... ";
		SR3 = SRExtract(setupfilenames[r].c_str(), varname);

		std::cout << endl << "Loading evaluation Sampler..." << endl;
		varname = "EvalSampProp";
		EvalSamp = SamplerExtract(setupfilenames[r].c_str(), varname, hlb);
		std::cout << endl << "Loading observable Operators..." << endl;
		for (int o = 0; o < OpNames.size(); o++)
		{
			EvalOps.push_back(OperatorExtract(setupfilenames[r].c_str(), OpNames[o]));
		}
		EvalVals.resize(EvalOps.size());

		std::cout << endl << "Loading Ansatz object..." << endl;
		varname = "AnsProp";
		Ans = AnsatzExtract(wfuncfilenames[r].c_str(), varname);
		varname = "StartParams";
		VectorMatSave(Ans->ParamList(), wfuncfilenames[r].c_str(), varname);

		// std::cout << "Testing common functions... " << endl;
		// Config* cfg = hlb->RandomCfg();
		// Ans->PrepPsi(cfg);
		// Diff diff = hlb->PropMove(cfg);
		// vector<int> cfg_vec = cfg->FullCfg();
		// VectorXd config = VectorXd::Zero(cfg_vec.size());
		// std::cout << "Configuration: ";
		// for (int c = 0; c < cfg_vec.size(); c++)
		// {
		// 	config(c) = cfg_vec[c];
		// 	std::cout << cfg_vec[c] << " ";
		// }
		// ParamSave_ASCII(config, "configuration.txt");
		// std::cout << endl << "Difference: " << endl;
		// for (int d = 0; d < diff.num; d++)
		// {
		// 	std::cout << "Position " << (d + 1) << ": " << diff.pos[d] << ", value " << diff.val[d] << endl;
		// }
		// std::cout << "Ratio: " << Ans->PsiRatio(diff) << endl;

		// VectorXd dLogp = Ans->LogDeriv(cfg);
		// std::cout << "Derivative vector: ";
		// for (int d = 0; d < dLogp.size(); d++)
		// {
		// 	std::cout << dLogp(d) << " ";
		// }
		// std::cout << endl;
		// ParamSave_ASCII(dLogp, "logderiv.txt");

		tstart = high_resolution_clock::now();
		EvalSamp->EvalSample(Ans, EneGS, EvalOps, EvalVals);
		tfinish = high_resolution_clock::now();
		trun = duration_cast<duration<double>>(tfinish - tstart);
		EvalTime = trun.count();
		std::cout << endl << "Evaluation sampling completed in " << trun.count() << " seconds." << endl;
		std::cout << "Energy evaluated by Monte Carlo sampling: " << EneGS << endl;
		varname = "EvalTime";
		ScalarMatSave(EvalTime, wfuncfilenames[r].c_str(), varname);
		varname = "EneGS";
		ScalarMatSave(EneGS, wfuncfilenames[r].c_str(), varname);
		for (int o = 0; o < ValNames.size(); o++)
		{
			MatrixMatSave(EvalVals[o], wfuncfilenames[r].c_str(), ValNames[o]);
		}

		tstart = high_resolution_clock::now();
		std::cout << endl << "Performing optimisation passes with SRObj1..." << endl;
		std::cout << "Number of samples per step: " << Samp1->NumSamp() << endl;
		EnIter = SR1->Optimise(Ans, Samp1);
		varname = "ParamPass1";
		VectorMatSave(Ans->ParamList(), wfuncfilenames[r].c_str(), varname);
		varname = "EnIter1";
		VectorMatSave(EnIter, wfuncfilenames[r].c_str(), varname);
		std::cout << endl << "Performing optimisation passes with SRObj2..." << endl;
		std::cout << "Number of samples per step: " << Samp2->NumSamp() << endl;
		EnIter = SR2->Optimise(Ans, Samp2);
		varname = "ParamPass2";
		VectorMatSave(Ans->ParamList(), wfuncfilenames[r].c_str(), varname);
		varname = "EnIter2";
		VectorMatSave(EnIter, wfuncfilenames[r].c_str(), varname);
		std::cout << endl << "Performing optimisation passes with SRObj3..." << endl;
		std::cout << "Number of samples per step: " << Samp3->NumSamp() << endl;
		EnIter = SR3->Optimise(Ans, Samp3);
		varname = "ParamPass3";
		VectorMatSave(Ans->ParamList(), wfuncfilenames[r].c_str(), varname);
		varname = "EnIter3";
		VectorMatSave(EnIter, wfuncfilenames[r].c_str(), varname);
		std::cout << endl << "Optimisation of " << wfuncfilenames[r] << " with " << setupfilenames[r] << " complete." << endl << endl;
		tfinish = high_resolution_clock::now();
		trun = duration_cast<duration<double>>(tfinish - tstart);
		RunTime = trun.count();
		varname = "RunTime";
		ScalarMatSave(RunTime, wfuncfilenames[r].c_str(), varname);
	}
	return 0;
}