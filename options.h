#ifndef OPTIONS_H
#define OPTIONS_H

#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>

using namespace std;

class options
{
public:
	int nfiles;
	int stride;
	vector<string> variables;
	string dataset;

	options(int argc, char* argv[])
	{
		try
		{
			// ************* Reader and Description ************* //
			TCLAP::CmdLine inp("Dymode, copyrighted for money", ' ', "0.1a");

			// Number of .h5 files to read
			TCLAP::ValueArg<int> nfilesArg("n", "nfiles", "Number of files to read", false /*req*/, 1/*default*/, "int", inp);

			// Stride between snapshots
			TCLAP::ValueArg<int> nskipstepArg("s", "stride", "Step between snapshots to read (read every other s snapshots", false /*req*/, 1/*default*/, "int", inp);

			// Dataset name within the h5 files to read from
			TCLAP::ValueArg<string> datasetnameArg("d", "dataset", "dataset name within the HDF file(s)", false /*req*/, "snapshots_T"/*default*/, "string", inp);

			// Name of H5 files
			TCLAP::ValueArg<string> filenameArg("f", "filename", "name of the data-file(s), without trailing number (rootname)", false /*req*/, "D:/DMD/DMD/x64/NNDEB/Re350_oscillating"/*default*/, "string", inp);
			
			// Variables to keep
			TCLAP::ValueArg<string> variablesArg("i", "variables", "name of the input variable(s) to keep in the snapshot matrix before starting the DMD. A name must be provided for each variable present in the disk data, separated by commas. Use 'null' in order to not use a variable. For example, if the data on disk contains the variables u, v, w, p but you only want to use u and p, use --variables u,null,null,p", false /*req*/, "null"/*default*/, "string", inp);


			// ************* Parse the argv array ************* //
			inp.parse(argc, argv);

			// ************* Store the arguments ************** //
			nfiles = nfilesArg.getValue();
			stride = nskipstepArg.getValue();

			boost::split(variables, variablesArg.getValue(), boost::is_any_of(","));
			dataset = datasetnameArg.getValue();
		}
		catch (TCLAP::ArgException &e)  // catch any exceptions
		{
			std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		}
	}
};

#endif //OPTIONS_H