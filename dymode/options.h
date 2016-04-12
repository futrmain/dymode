#ifndef OPTIONS_H
#define OPTIONS_H

//#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include "yaml-cpp/yaml.h"

using namespace std;

class options
{
public:
	int nfiles;
	int stride;
	vector<string> variables;
	string dataset;
	string filename;
	string outdir;
	bool dispResiduals;
	peigen::EigenMethod eigSolver;
	int nmodes;
	int npod;
	string geofile;
	int nsingulars;
	int sblock;
	sort_method sortMeth;
	string result;

	options(int argc, char* argv[])
	{
		try
		{
//cout << "I'm good" << endl;
			if (argc < 2)
			{
//cout << "Oopsy" << endl;
				result = "The path to YAML file is missing.";
			}
			else
			{
//cout << "Srious business" << endl;
				//cout << argv[1] << endl;
				//cout << "coutchi" << endl;
//cout << "coutchi 2" << endl;
				YAML::Node config = YAML::LoadFile(argv[1]);
				//cout << "fuck you" << endl;
				//cout << config["input"]["filename"].as<string>() << endl;

				YAML::Node subnode = config["input"];
//cout << subnode.size() << endl;
//cout << "OK" << endl;
/*for (YAML::const_iterator it=subnode.begin();it!=subnode.end();++it)
{
	//cout << i << endl;
	cout << it->first.as<std::string>() << endl;
	cout << it->second.as<std::string>() << endl;
}*/
				filename = subnode["filename"].as<std::string>();
				//cout << filename << endl;

				YAML::Node dsnode = subnode["datasets"];
//cout << "number of datasets: " << dsnode.size() << endl;					
				if (subnode["datasets"].size() < 1)	
				{
//cout << "You must specify at least one dataset." << endl;
					result = "You must specify at least one dataset.";
				}
				else
				{
//cout << "number of datasets: " << subnode["datasets"].size() << endl;
					for (int i = 0; i < subnode["datasets"].size(); ++i)
					{
//cout << subnode["datasets"][i].as<std::string>() << endl;
						variables.push_back(subnode["datasets"][i].as<std::string>());
					}
				}

				if (subnode["nfiles"])
				{
					nfiles = subnode["nfiles"].as<int>();
				}
				else
				{
					nfiles = 1;
				}
				//cout << nfiles << endl;

				if (subnode["stride"])
				{
					stride = subnode["stride"].as<int>();
				}
				else
				{
					stride = 1;
				}
				//cout << stride << endl;

				if (config["computation"])
				{
					subnode = config["computation"];

					if (subnode["block"])
					{
						sblock = subnode["block"].as<int>();
					}
					else
					{
						sblock = 6;
					}

					if (subnode["eigen"])
					{
						if (subnode["eigen"].as<string>() == "EigSerial")
							eigSolver = EigSerial;
						else if (subnode["eigen"].as<string>() == "EigHess")
							eigSolver = EigHess;
						else if (subnode["eigen"].as<string>() == "EigSchur")
							eigSolver = EigSchur;
						else // unrecognized option
							eigSolver = EigHess;
					}
					else
					{
						eigSolver = EigHess;
					}

					if (subnode["residuals"])
					{
						dispResiduals = subnode["residuals"].as<bool>();
					}
					else
					{
						dispResiduals = false;
					}

					if (subnode["singulars"])
					{
						nsingulars = subnode["singulars"].as<int>();
					}
					else
					{
						nsingulars = 0;
					}
				}
				else
				{
					sblock = 6;
					eigSolver = EigHess;
					dispResiduals = false;
					nsingulars = 0;
				}

				subnode = config["output"];
//cout << "tada" << endl;
				outdir = subnode["outdir"].as<string>();
//cout << outdir << endl;
				if (subnode["sorting"])
				{
					sortMeth.set(subnode["sorting"].as<string>());
				}
				else
				{
					sortMeth.set("energy,median");
				}

				if (subnode["pod"])
				{
					npod = subnode["pod"].as<int>();
				}
				else
				{
					npod = 0;
				}

				if (subnode["dmd"])
				{
					nmodes = subnode["dmd"].as<int>();
				}
				else
				{
					nmodes = 0;
				}	
//cout << "/else" << endl;			
			}
		}
		catch(YAML::ParserException& e)
		{
			std::cout << "erreur!\n";
			std::cout << e.what() << "\n";
		}
//cout << "big boom" << endl;
	}
};


		

#endif //OPTIONS_H
