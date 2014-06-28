#ifndef GEOINDMD_H
#define GEOINDMD_H

#include "mpi.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include "Eigen/Dense"
#include <fstream>


using namespace std;
using namespace Eigen;

struct gold_part
{
	int number;
	vector<int> nelements;
	vector<string> telements;
};

class geofilereader
{
private:
	const map<string, int> np = { { "point", 1 }, { "gpoint", 1 }, { "bar2", 2 }, { "g_bar2", 2 }, { "bar3", 3 }, { "g_bar3", 3 }, { "tria3", 3 }, { "g_tria3", 3 }, { "tria6", 6 }, { "g_tria6", 6 }, { "quad4", 4 }, { "g_quad4", 4 }, { "quad8", 8 }, { "g_quad8", 8 }, { "tetra4", 4 }, { "g_tetra4", 4 }, { "tetra10", 10 }, { "g_tetra10", 10 }, { "pyramid5", 5 }, { "g_pyramid5", 5 }, { "pyramid13", 13 }, { "g_pyramid13", 13 }, { "penta6", 6 }, { "g_penta6", 6 }, { "penta15", 15 }, { "g_penta15", 15 }, { "hexa8", 8 }, { "g_hexa8", 8 }, { "hexa20", 20 }, { "g_hexa20", 20 }, { "nsided", -1 }, { "g_nsided", -1 }, { "nfaced", -1 }, { "g_nfaced", -1 } };

public:
	string filename;
	vector<gold_part> parts;
	

	geofilereader() {};
	~geofilereader() {};

	geofilereader(const string& geofile) : filename(geofile)
	{
		parse(geofile);
	};

	void parse(const string& geofile)
	{
		char line[81], subLine[81], nameline[81];
		line[80] = '\0'; subLine[80] = '\0'; nameline[80] = '\0';
#ifdef _WIN32
		ifstream input = ifstream(geofile.c_str(), ios::in | ios::binary);
#else
		ifstream input = ifstream(geofile.c_str(), ios::in);
#endif

		int lineRead;

		//FIXME should switch back to no-geo file behavior instead of crashing, especially if it occurs after expensive computation...
		assert(!input.fail() && "COULD NOT OPEN THE GEOMETRY FILE PROVIDED");

		// Skip the first 3 lines
		input.read(line, 80);
		input.read(line, 80);
		input.read(line, 80);
		cout << BLACS::myrank << " : " << line << endl;

		bool NodeIdsListed, ElementIdsListed;
		// Read the node id and element id lines.
		input.read(line, 80);
		cout << BLACS::myrank << " : " << line << endl;
		std::sscanf(line, " %*s %*s %s", subLine);
		cout << BLACS::myrank << " : " << subLine << endl;
		if (strncmp(subLine, "given", 5) == 0)
		{
			NodeIdsListed = true;
		}
		else if (strncmp(subLine, "ignore", 6) == 0)
		{
			NodeIdsListed = true;
		}
		else
		{
			NodeIdsListed = false;
		}

		input.read(line, 80);
		cout << BLACS::myrank << " : " << line << endl;
		std::sscanf(line, " %*s %*s %s", subLine);
		cout << BLACS::myrank << " : " << subLine << endl;
		if (strncmp(subLine, "given", 5) == 0)
		{
			ElementIdsListed = true;
		}
		else if (strncmp(subLine, "ignore", 6) == 0)
		{
			ElementIdsListed = true;
		}
		else
		{
			ElementIdsListed = false;
		}

		input.read(line, 80); // "extents" or "part"
		cout << BLACS::myrank << " : " << line << endl;
		if (strncmp(line, "extents", 7) == 0)
		{
			// Skipping the extents.
			//input.seekg(6 * sizeof(float), ios::cur);
			input.ignore(6 * sizeof(float)); // "part"
			input.read(line, 80); // "part"
			cout << BLACS::myrank << " : " << line << endl;
		}

		while (input.good() && strncmp(line, "part", 4) == 0)
		{
			gold_part curr_part;
			input.read((char*)&(curr_part.number), sizeof(int));

			assert(curr_part.number >= 0 && curr_part.number < 65536 && "The part numer found in the geometry file is not correct.");

			input.read(line, 80); // part description line
			cout << BLACS::myrank << " : " << line << endl;

			input.read(line, 80); // coordinates
			cout << BLACS::myrank << " : " << line << endl;
			assert(strncmp(line, "block", 5) && "Block syntax is not supported.");
			assert(!strncmp(line, "coordinates", 11) && "Expected to find 'coordinates' in geo file");

			int nel, nnodes;
			string tel;

			input.read((char*)(&nnodes), sizeof(int)); // Number of nodes
			cout << BLACS::myrank << "nnodes : " << nnodes << endl;

			if (NodeIdsListed)
			{	// Skip element IDs
				input.ignore(nnodes * sizeof(int));
			}
			input.ignore(3 * nnodes * sizeof(float)); // Skip node coordinates

			input.read(line, 80); // element type
			cout << BLACS::myrank << " : " << line << endl;
			std::sscanf(line, " %s", subLine);
			auto it = np.find(subLine);
			cout << BLACS::myrank << " : " << subLine << endl;
			while (input.good() && it != np.end())
			{
				int npoints; 
				npoints = it->second;

				curr_part.telements.push_back(string(subLine));

				input.read((char*)(&nel), sizeof(int)); // Number of elements
				curr_part.nelements.push_back(nel);

				

				if (ElementIdsListed)
				{	// Skip element IDs
					input.ignore(nel * sizeof(int));
				}
				if (npoints > 0)
				{
					input.ignore(npoints * nel * sizeof(int)); // skip elements
				}
				else // n-sided hell
				{
					// read the number of points for each element
					int* nnp = new int[nel];
					input.read((char*)nnp, nel * sizeof(int));

					int total = 0;
					for (int k = 0; k < nel; ++k)
					{
						total += nnp[k];
					}



					input.ignore(total * sizeof(int)); // skip elements
					delete[] nnp;
				}

				input.read(line, 80); // element type or "part"
				cout << BLACS::myrank << " : " << line << endl;
				std::sscanf(line, " %s", subLine);
				it = np.find(subLine);
			}

			parts.push_back(curr_part);
			//input.read(line, 80); // "part"
			cout << BLACS::myrank << " : " << line << endl;
		}

		input.close();
	};
};

#endif // GEOINDMD_H