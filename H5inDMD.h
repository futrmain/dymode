#ifndef H5INDMD_H
#define H5INDMD_H

#include "hdf5.h"
#include "SharedMatrix.h"
#include "mpi.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>

namespace phdfp
{
	using namespace std;
	using namespace peigen;

	struct hyperslab2D
	{
		hsize_t offset[2];
		hsize_t stride[2];
		hsize_t count[2];
		hsize_t block[2];
	};

	class datasetreader
	{
	public:
		Matrix<double, Dynamic, Dynamic> input_data;
		
		hsize_t virtual_dims[2];
		int nfiles;	// NOTE!!nfiles should be places before firstname because it needs to have been initialized when firstname is initializing 
		string rootname;
		string firstname;

		Matrix<Matrix<int, 1, Dynamic>, Dynamic, 1> snapshots_per_file;

		// Constructor
		datasetreader(string fname) : firstname(fname), nfiles(1) 
		{}

		datasetreader(int num_files, string root_name) : nfiles(num_files), rootname(root_name), firstname(filename(1))
		{
			//cout << BLACS::myrank << " will resize snapshots_per_file " << BLACS::numproc << "x" << 1 << "" << endl;
			snapshots_per_file.resize(BLACS::numproc, 1);
		}

		// Destructor
		~datasetreader() { /*H5Sclose(virtual_filespace);*/ }

		// Members
		string filename(int i)
		{
			//std::cout << "hello" << std::endl;
			ostringstream  file;
			file << rootname << std::setw(4) << setfill('0') << i + 0 << ".h5";
			//std::cout << (nfiles > 1 ? file.str() : file.str()) << std::endl;
			return (nfiles > 1 ? file.str() : file.str());
		}

		void read(string dataset_name);
		void getextents(string fname, string dataset_name, hsize_t * dims);

		void createShared(int rblock, int cblock);
	};


	void datasetreader::getextents(string fname, string dataset_name, hsize_t * dims)
	{
		hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);
		hid_t filespace = H5Dget_space(dataset);

		assert((H5Sget_simple_extent_dims(filespace, dims, NULL) == 2) && "THE RANK OF THE DATASET IS NOT 2");
		H5Sget_simple_extent_dims(filespace, dims, NULL);

		H5Sclose(filespace);
		H5Dclose(dataset);
		H5Fclose(file);
	}

	void datasetreader::read(string dataset_name)
	{
		// Step 0: Figure out how many files for each process
		for (int proc = 0; proc < nfiles % BLACS::numproc; ++proc)
		{
			//cout << BLACS::myrank << " ahah " << proc << " will resize " << 1 << "x" << floor(nfiles / BLACS::numproc) + 1 << "" << endl;
			snapshots_per_file(proc).resize(1, floor(nfiles / BLACS::numproc) + 1); 
		}
		for (int proc = nfiles % BLACS::numproc; proc < BLACS::numproc; ++proc)
		{
			//cout << BLACS::myrank << " and " << proc << " will resize " << 1 << "x" << floor(nfiles / BLACS::numproc) << "" << endl;
			snapshots_per_file(proc).resize(1, floor(nfiles / BLACS::numproc));
		}

		// Step 1: Figure out how much data will be read by this process
		int dimension_r = 0;
		int dimension_c = 0;

		
		for (int k = 0; k < snapshots_per_file(BLACS::myrank).size(); ++k)
		{
			int filenum = BLACS::myrank + 1 + k * BLACS::numproc;
			string fname = filename(filenum);
			
			hsize_t dims[2];
			getextents(fname, dataset_name, dims);

			assert(((dimension_r == 0) || (dimension_r == dims[1])) && "NUMBER OF ROWS MISMATCH FROM ONE FILE TO ANOTHER");
			dimension_r = dims[1];
			snapshots_per_file(BLACS::myrank)(k) = dims[0];
		}

		// Step 2: Allocate the memory
		input_data.resize(dimension_r, snapshots_per_file(BLACS::myrank, 0).sum());

		cout << BLACS::myrank << " will read " << input_data.rows() << "x" << input_data.cols() << "data" << endl;

		// Step 3: Read
		size_t dims[2];
		dims[0] = input_data.rows();
		dims[1] = input_data.cols();
		hsize_t virtual_filespace = H5Screate_simple(/*rank*/ 2, dims, NULL);

		hyperslab2D sel;
		sel.offset[0] = 0;
		sel.offset[1] = 0;
		sel.stride[0] = 1;
		sel.stride[1] = 1;
		sel.count[0] = 1;
		sel.count[1] = 1;
		for (int k = 0; k < snapshots_per_file(BLACS::myrank).size(); ++k)
		{
			sel.block[0] = input_data.rows();
			sel.block[1] = snapshots_per_file(BLACS::myrank)(k);
			H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);

			int filenum = BLACS::myrank + 1 + k * BLACS::numproc;
			string fname = filename(filenum);
			//cout << fname << endl << flush;

			hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

			hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);

			herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, virtual_filespace, H5S_ALL, H5P_DEFAULT, input_data.data());

			H5Dclose(dataset);
			H5Fclose(file);
			sel.offset[1] += sel.block[1];
		}
		H5Sclose(virtual_filespace);

		// Step 4: Exchange how many snapshots in each file
		for (int proc = 0; proc < BLACS::numproc; ++proc)
		{
			MPI::COMM_WORLD.Bcast(snapshots_per_file(proc).data(), snapshots_per_file(proc).size(), MPI::DOUBLE, proc);
		}
		/*if (BLACS::myrank == 2)
		{
			for (int proc = 0; proc < BLACS::numproc; ++proc)
				cout << snapshots_per_file(proc) << ";" << endl << endl;
		}*/
			
	}


	void datasetreader::createShared(int rblock, int cblock)
	{
		// Step 1: Figure out the global index of all mmy columns
		MatrixXi global_indices(1,input_data.cols());
		int i_glob = 0;
		int i_loc = 0;

		//cout << ceil((float)nfiles / BLACS::numproc) << endl;
		for (int f = 0; f < floor(nfiles / BLACS::numproc); ++f)
		{
			for (int p = 0; p < BLACS::myrank; ++p)
			{
				i_glob += snapshots_per_file(p, 0)(0, f);
			}

			for (int c = 0; c < snapshots_per_file(BLACS::myrank, 0)(0, f); ++c)
			{
				global_indices(0, i_loc + c) = i_glob;
				++i_glob;
			}
			i_loc += snapshots_per_file(BLACS::myrank, 0)(0, f);

			for (int p = BLACS::myrank + 1; p < BLACS::numproc; ++p)
			{
				i_glob += snapshots_per_file(p, 0)(0, f);
			}
		}

		for (int proc = 0; proc < nfiles % BLACS::numproc; ++proc)
		{
			if (proc == BLACS::myrank)
			{
				for (int c = 0; c < snapshots_per_file(BLACS::myrank, 0)(0, ceil((float)nfiles / BLACS::numproc)-1); ++c)
				{
					global_indices(0, i_loc + c) = i_glob;
					++i_glob;
				}
			}
			i_glob += snapshots_per_file(proc, 0)(0, ceil((float)nfiles / BLACS::numproc)-1);
		}

			/*for (int c = 0; c < snapshots_per_file(BLACS::myrank,0)(0,f); ++c)
			{
				if (snapshots_per_file(BLACS::myrank, 0).size() > f)
				{
					global_indices(0, i_loc + c) = i_glob;
					++i_glob;
				}
			}
			i_loc += snapshots_per_file(BLACS::myrank,0)(0,f);

			for (int p = BLACS::myrank+1; p < BLACS::numproc; ++p)
			{
				if (snapshots_per_file(p,0).size() > f)
				{
					i_glob += snapshots_per_file(p,0)(0,f);
				}
			}*/
		

		cout << BLACS::myrank << " " << global_indices << endl << endl << endl;


	}


}	// end namespace phdfp
#endif // H5INDMD_H