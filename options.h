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
	string outdir;
	bool dispResiduals;
	peigen::EigenMethod eigSolver;
	int nmodes;
	string geofile;

	options(int argc, char* argv[])
	{
		try
		{
			// ************* Reader and Description ************* //
			TCLAP::CmdLine input("Dymode, copyrighted for money", ' ', "0.1a");

			// Number of .h5 files to read
			TCLAP::ValueArg<int> nfilesArg("n", "nfiles", "Number of files to read", false /*req*/, 1/*default*/, "int", input);

			// Stride between snapshots
			TCLAP::ValueArg<int> nskipstepArg("s", "stride", "Step between snapshots to read (read every other s snapshots", false /*req*/, 1/*default*/, "int", input);

			// Dataset name within the h5 files to read from
			TCLAP::ValueArg<string> datasetnameArg("d", "dataset", "dataset name within the HDF file(s)", false /*req*/, "snapshots_T"/*default*/, "string", input);

			// Name of H5 files
			TCLAP::ValueArg<string> filenameArg("f", "filename", "name of the data-file(s), without trailing number (rootname)", false /*req*/, "D:/DMD/DMD/x64/NNDEB/Re350_oscillating"/*default*/, "string", input);
			
			// Variables to keep
			TCLAP::ValueArg<string> variablesArg("i", "variables", "name of the input variable(s) to keep in the snapshot matrix before starting the DMD. A name must be provided for each variable present in the disk data, separated by commas. Use 'null' in order to not use a variable. For example, if the data on disk contains the variables u, v, w, p but you only want to use u and p, use --variables u,null,null,p", false /*req*/, "null"/*default*/, "string", input);

			// Output directory
			TCLAP::ValueArg<string> outdirArg("o", "outdir", "Output directory where result files are saved", false /*req*/, ""/*default*/, "string", input);

			// Display residuals switch
			TCLAP::SwitchArg residualsArg("r", "residuals", "If specified, residuals are computed after the key steps of the computation", input, false);

			// Method to use in ScaEigenSolver
			TCLAP::ValueArg<string> eigArg("e", "eigen", "String describing the level of parallelism to use when solving the eigen value problem. Possible values are \"EigSerial\": the whole problem is solved in serial; \"EigHess\": Hessenberg reduction is done in parallel, the rest is done in serial; \"EigSchur\": Everything is computed in parallel. This is bugged in MKL 11.", false /*req*/, "EigHess"/*default*/, "string", input);

			// Number of modes to print out
			TCLAP::ValueArg<int> nmodesArg("m", "modes", "Number of modes to save to disk", false /*req*/, 1/*default*/, "int", input);

			// Geometry file from Ensight Gold
			TCLAP::ValueArg<string> geoArg("g", "geo", "Geometry file from Ensight", false /*req*/, "D:/DMD/DMD/x64/NNDEB/dmd.geo"/*default*/, "string", input);

			// ************* Parse the argv array ************* //
			input.parse(argc, argv);

			// ************* Store the arguments ************** //
			nfiles = nfilesArg.getValue();
			stride = nskipstepArg.getValue();

			boost::split(variables, variablesArg.getValue(), boost::is_any_of(","));
			dataset = datasetnameArg.getValue();

			outdir = outdirArg.getValue();
			dispResiduals = residualsArg.getValue();

			if (eigArg.getValue() == "EigSerial")
				eigSolver = EigSerial;
			else if (eigArg.getValue() == "EigHess")
				eigSolver = EigHess;
			else if (eigArg.getValue() == "EigSchur")
				eigSolver = EigSchur;
			else // unrecognized option
				eigSolver = EigHess;

			nmodes = nmodesArg.getValue();

			geofile = geoArg.getValue();
		}
		catch (TCLAP::ArgException &e)  // catch any exceptions
		{
			std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		}
	}
};

#endif //OPTIONS_H