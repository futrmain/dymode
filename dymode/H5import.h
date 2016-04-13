#ifndef H5IMPORT_H
#define H5IMPORT_H

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

	struct range1D
	{
		hsize_t beg;
		hsize_t end;
	};

	class datasetreader
	{
	public:
		Matrix<double, Dynamic, Dynamic, ColMajor> input_data;
		
		int virtual_dims[2];
		int nfiles;	// NOTE!!nfiles should be placed before firstname because it needs to have been initialized when firstname is initializing 
		string rootname;
		string firstname;
		int Np; // Number of points in the geometry (points per variable)
		int nvars; // Number of variables to include in the DMD
		map<string, range1D> ranges;

		Matrix<Matrix<int, 1, Dynamic>, Dynamic, 1> snapshots_per_process;

		// Constructor
		datasetreader(string fname) : nfiles(1), firstname(fname)
		{}

		datasetreader(int num_files, string root_name) : nfiles(num_files), rootname(root_name), firstname(filename(1))
		{
			//cout << BLACS::myrank << " will resize snapshots_per_file " << BLACS::numproc << "x" << 1 << "" << endl;
			snapshots_per_process.resize(BLACS::numproc, 1);
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

		//void read(string dataset_name);
		void read(vector<string> datasets);
		void getextents(string fname, string dataset_name, hsize_t * dims);

		peigen::SharedMatrix<MatrixXd> createShared(int rblock, int cblock, int nskip = 1);
	};


	inline void datasetreader::getextents(string fname, string dataset_name, hsize_t * dims)
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

	inline void datasetreader::read(vector<string> variables)
	{
		// Step 0: Figure out how many files for each process
		int my_num_files = floor(nfiles / BLACS::numproc);

		if (BLACS::myrank < nfiles % BLACS::numproc)
		{
			++my_num_files;
		}


		/*
		// Step 1.1: Figure out how many variables will effectively be read
		nvars = 0;
		for (string v : variables)
		{
			if (!(v == "null"))
			{
				++nvars;
			}
		}
		*/

		int dimension_r = 0;
		int dimension_c = 0;

		for (int field = 0; field < variables.size(); ++field)
		{
			ranges[variables[field]].beg = 0;
			ranges[variables[field]].end = 0;
		}

		MatrixXi snapshots_per_file(my_num_files, 1);

		// Step 1.2: Figure out how much data will be read by this process
		for (int f = 0; f < my_num_files; ++f)
		{
			int filenum = BLACS::myrank * floor(nfiles / BLACS::numproc) + min(BLACS::myrank, nfiles % BLACS::numproc) + f + 1;
			string fname = filename(filenum);

			int idx = 0;
			dimension_r = 0;
			int ncols_per_field = 0;

			for (int field = 0; field < variables.size(); ++field)
			{
				hsize_t dims[2];
				getextents(fname, variables[field], dims);

				int rbeg, rend;
				rbeg = idx;
				rend = idx + dims[1] - 1;
				idx += dims[1];
		
				assert(((ranges[variables[field]].beg == 0) || (ranges[variables[field]].beg == rbeg)) && "NUMBER OF ROWS MISMATCH FROM ONE FILE TO ANOTHER.");	
				assert(((ranges[variables[field]].end == 0) || (ranges[variables[field]].end == rend)) && "NUMBER OF ROWS MISMATCH FROM ONE FILE TO ANOTHER.");	

				assert(((ncols_per_field == 0) || (ncols_per_field == dims[0])) && "NUMBER OF SNAPSHOTS MISMATCH FROM ONE DATASET TO ANOTHER.");	

				ranges[variables[field]].beg = rbeg;
				ranges[variables[field]].end = rend;

				snapshots_per_file(f, 0) = dims[0];
				ncols_per_field = dims[0];
				dimension_r += dims[1];
			}
			dimension_c += ncols_per_field;
			//cout << BLACS::myrank << " reads " << dimension_c << endl;
		}
		
		//cout << BLACS::myrank << " Np =  " << Np << "... " << "dimension_r = " << dimension_r << "variables.size()" << variables.size() << endl;
		virtual_dims[0] = dimension_r;


		// Step 2: Allocate the memory
		snapshots_per_process(BLACS::myrank, 0).resize(1, dimension_c);
		input_data.resize(dimension_r, dimension_c);

		if (input_data.rows() * input_data.cols() > 0) // Keep quiet if no data to read
		{
			//cout.precision(6);
			cout << "(" << BLACS::myrank << ") will read " << input_data.rows() << "x" << input_data.cols() << " values\t(" << input_data.rows() * input_data.cols() * 8 * 0.000001 << " MiB)" << endl;
		}

		//// Step 3: Read
		/*Note: From here the dimensions are switched to account for the fact that data is stored transposed on disk*/
		hsize_t dims[2];
		dims[1] = input_data.rows();
		dims[0] = input_data.cols();
		hid_t virtual_filespace = H5Screate_simple(/*rank*/ 2, dims, NULL);

		hyperslab2D sel;
		sel.stride[0] = 1;
		sel.stride[1] = 1;
		sel.count[0] = 1;
		sel.count[1] = 1;
		sel.offset[0] = 0;
		for (int f = 0; f < my_num_files; ++f)
		{
			for (int field = 0; field < variables.size(); ++field)
			{
				sel.offset[1] = ranges[variables[field]].beg;
				sel.block[1] = ranges[variables[field]].end - ranges[variables[field]].beg + 1; //input_data.rows();
				sel.block[0] = snapshots_per_file(f, 0);
				H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);
			
cout << "selected a block in virtual space " << sel.offset[0]<<","<<sel.offset[1] <<", " << sel.block[0]<<","<<sel.block[1] <<endl;
				int filenum = BLACS::myrank * floor(nfiles / BLACS::numproc) + min(BLACS::myrank, nfiles % BLACS::numproc) + f + 1;	
				string fname = filename(filenum);
				//cout << fname << endl << flush;
	
				hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	
				hid_t dataset = H5Dopen2(file, variables[field].c_str(), H5P_DEFAULT);
				hsize_t fdim[2];
				fdim[0] = snapshots_per_file(f, 0);
				fdim[1] = ranges[variables[field]].end - ranges[variables[field]].beg + 1;
				hsize_t file_filespace = H5Screate_simple(/*rank*/ 2, fdim, NULL);
				hyperslab2D fsel; 
				fsel.block[1] = ranges[variables[field]].end - ranges[variables[field]].beg + 1;
				fsel.block[0] = snapshots_per_file(f, 0);
				fsel.offset[0] = 0;
				fsel.offset[1] = 0;
				fsel.stride[0] = 1;
				fsel.stride[1] = 1;
				fsel.count[0] = 1;
				fsel.count[1] = 1;
				auto op = H5S_SELECT_SET;
				
				H5Sselect_hyperslab(file_filespace, op, fsel.offset, fsel.stride, fsel.count, fsel.block);
//cout << "selected a block in file space " << fsel.offset[0]<<","<<fsel.offset[1] <<", " << fsel.block[0]<<","<<fsel.block[1] <<endl;
					
				herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, virtual_filespace, file_filespace, H5P_DEFAULT, input_data.data());
	
				H5Dclose(dataset);
				H5Fclose(file);				
			}
			sel.offset[0] += sel.block[0];
		}
		H5Sclose(virtual_filespace);
//cout << "step 4 beg " <<endl;
		// Step 4: Exchange how many snapshots for each process
		for (int proc = 0; proc < BLACS::numproc; ++proc)
		{
			int n_snaps = snapshots_per_process(proc, 0).size();
			MPI::COMM_WORLD.Bcast(&n_snaps, 1, MPI::INT, proc);
			snapshots_per_process(proc, 0).resize(n_snaps);
		}
//cout << "step 4 end" <<endl;
		// Step 5: Broadcast the number of rows for processes that don't have any file to open
		// It is assumed that at least proc 0 has a file to read...
		MPI::COMM_WORLD.Bcast(virtual_dims, 1, MPI::INT, 0);

		Np = virtual_dims[0];
//cout << Np <<endl;
		/*if (BLACS::myrank == 2)
		{
		for (int proc = 0; proc < BLACS::numproc; ++proc)
		cout << snapshots_per_file(proc) << ";" << endl << endl;
		}*/

	}

	
	inline peigen::SharedMatrix<MatrixXd> datasetreader::createShared(int rblock, int cblock, int nskip)
	{
		// Step 1: Figure out the global index of all my snapshots
		int i_glob = 0;

		for (int p = 0; p < BLACS::numproc; ++p)
		{
			for (int c = 0; c < snapshots_per_process(p, 0).size(); ++c)
			{
				// If the column should be read, tag it with a virtual global index
				if (i_glob % nskip == 0)
				{
					snapshots_per_process(p, 0)(0, c) = floor(i_glob / nskip);
				}
				// If it shouldn't be in the SharedMatrix, tag it with -1
				else
				{
					snapshots_per_process(p, 0)(0, c) = -1;
				}
				i_glob++;
			}
			//cout << BLACS::myrank << " showing for proc " << p << ": " << snapshots_per_process(p, 0) << endl << endl << endl;
		}


		// Step 2: Prepare send buffer
		Matrix<Matrix<double, Dynamic, Dynamic>, Dynamic, Dynamic> SendBuff(BLACS::grid_rows, BLACS::grid_cols);
		

		// Step 2.1: Figure out the size of each buffer
		MatrixXi cols2send = MatrixXi::Zero(1, BLACS::grid_cols);
		for (int c = 0; c < snapshots_per_process(BLACS::myrank, 0).size(); ++c)
		{
			int global_index = snapshots_per_process(BLACS::myrank, 0)(0, c);
			if (global_index >= 0)
			{
				cols2send(0, BLACS::indxg2p(global_index, cblock, BLACS::grid_cols))++;
			}
		}
		//cout << BLACS::myrank << " cols2send " << cols2send << endl;

		//MatrixXi rows2send = MatrixXi::Zero(BLACS::grid_rows, 1);
		Matrix<Matrix<int, 2, 1>, Dynamic, Dynamic> offsets(BLACS::grid_rows, BLACS::grid_cols);
		for (int r = 0; r < BLACS::grid_rows; ++r)
		{
			int rows2send = BLACS::numroc_(virtual_dims, &rblock, &r, BLACS::iZERO, &(BLACS::grid_rows));
			for (int pcol = 0; pcol < BLACS::grid_cols; ++pcol)
			{
				SendBuff(r, pcol).resize(rows2send, cols2send(0, pcol));
				offsets(r, pcol) = Matrix<int, 2, 1>::Zero();
			}
		}


		// Step 2.2: Fill the Send buffers

		for (int c = 0; c < snapshots_per_process(BLACS::myrank, 0).size(); ++c)
		{
			int global_index = snapshots_per_process(BLACS::myrank, 0)(0, c);
			if (global_index >= 0)
			{
				int pcol = BLACS::indxg2p(global_index, cblock, BLACS::grid_cols);
				int n_full_blocks = floor(virtual_dims[0] / rblock);

				for (int block = 0; block < n_full_blocks; ++block)
				{
					SendBuff(block % BLACS::grid_rows, pcol).block(offsets(block % BLACS::grid_rows, pcol)(0, 0), offsets(block % BLACS::grid_rows, pcol)(1, 0), rblock, 1)
						= input_data.block(block*rblock, c, rblock, 1);

					offsets(block % BLACS::grid_rows, pcol)(0, 0) += rblock;
				}

				SendBuff(n_full_blocks % BLACS::grid_rows, pcol).block(offsets(n_full_blocks % BLACS::grid_rows, pcol)(0, 0), offsets(n_full_blocks % BLACS::grid_rows, pcol)(1, 0), virtual_dims[0] % rblock, 1)
					= input_data.block(n_full_blocks*rblock, c, virtual_dims[0] % rblock, 1);

				// Next column for offsets
				for (int r = 0; r < BLACS::grid_rows; ++r)
				{
					offsets(r, pcol)(0, 0) = 0;
					offsets(r, pcol)(1, 0)++;
				}
			}
		}

		
		
		

		// Step 2.3: Clear source data to save space
		input_data.resize(0, 0);

		/*
		for (int rank = 0; rank < BLACS::numproc; ++rank)
		{
			if (BLACS::myrank == rank)
			{
				//cout << "The send buffers on " << peigen::BLACS::myrank << " are " << endl;
				for (int prow = 0; prow < BLACS::grid_rows; ++prow)
				{
					for (int pcol = 0; pcol < BLACS::grid_cols; ++pcol)
					{
						cout << prow << "x" << pcol << endl << SendBuff(prow, pcol) << endl << endl;
					}
				}
			}
			MPI::COMM_WORLD.Barrier();
		}
		*/
		
		// Step 3: Create the SharedMatrix, will be the receive buffer
		peigen::SharedMatrix<MatrixXd> S;
		if (BLACS::active)
		{
			S.resize(virtual_dims[0], ceil((float)i_glob / nskip), rblock, cblock);
		}
		else
		{
			S.local_matrix.resize(1, 1);
			// Processes not on the grid keep an uninitialized SharedMatrix
		}
		//S.printDetails();
		//cout << BLACS::myrank << " i_glob " << i_glob << endl;

		//cout << "The shared Matrix BEFORE on " << peigen::BLACS::myrank << " is " << endl << S.local_matrix << endl;

		// Step 4: Prepare to receive into the Shared Matrix
		int offset = 0;
		for (int p = 0; p < BLACS::numproc; ++p)
		{
			//How many snapshots will I receive from p?
			int nsnaps_from_p = 0;
			for (int s = 0; s < snapshots_per_process(p, 0).size(); ++s)
			{
				if (snapshots_per_process(p, 0)(0, s) >= 0)
				{
					int virtualGlobalIndex = snapshots_per_process(p, 0)(0, s);
					if (BLACS::indxg2p(virtualGlobalIndex, cblock, BLACS::grid_cols) == BLACS::mycol)
					{
						nsnaps_from_p += 1;
					}
				}
			}
			//cout << BLACS::myrank << " receive from " << p << " this many columns: " << nsnaps_from_p << endl;
			MPI::COMM_WORLD.Irecv(S.localData() + offset * S.local_matrix.rows(), S.local_matrix.rows() * nsnaps_from_p, MPI::DOUBLE, p, p);
			offset += nsnaps_from_p;
		}
		
		MPI::COMM_WORLD.Barrier();
		// Step 5: Send the data to everyone
		for (int prow = 0; prow < BLACS::grid_rows; ++prow)
		{
			for (int pcol = 0; pcol < BLACS::grid_cols; ++pcol)
			{
				int dest = BLACS::peigen_pnum(prow, pcol);
				MPI::COMM_WORLD.Send(SendBuff(prow, pcol).data(), SendBuff(prow, pcol).rows() * SendBuff(prow, pcol).cols(), MPI::DOUBLE, dest, BLACS::myrank);
			}
		}
		MPI::COMM_WORLD.Barrier();
		//(const void* buf, int count, const Datatype& datatype, int dest, int tag) const; 
		//BLACS::COMM_ACTIVE.Recv(localData(), local_matrix.size(), MPIType(), source, 1);

		//S.printDetails();
		return S;
	}
	

}	// end namespace phdfp
#endif // H5IMPORT_H
