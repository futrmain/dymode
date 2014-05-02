#ifndef PHDFP_DATASET
#define PHDFP_DATASET

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
		hsize_t dims[2];
		hsize_t virtual_dims[2];
		int nfiles;	// NOTE!!nfiles should be places before firstname because it needs to have been initialized when firstname is initializing 
		string rootname;
		string firstname;
		int rblock, cblock;

		

		SharedMatrix<MatrixXd> datamat;
		vector<hyperslab2D> selection;
		hid_t virtual_filespace;


		// Constructor
		datasetreader(string fname, int rowblock = peigen::BLACS::rblock, int colblock = peigen::BLACS::cblock) : firstname(fname), rblock(rowblock), cblock(colblock), nfiles(1) {}

		datasetreader(int num_files, string root_name, int rowblock = peigen::BLACS::rblock, int colblock = peigen::BLACS::cblock) : nfiles(num_files), rootname(root_name), rblock(rowblock), cblock(colblock), firstname(filename(1)) 
		{

			/*std::cout << "Number of files: " << num_files << std::endl 
				<< "Rootname: " << rootname << std::endl
				<< "RBLOCK: " << rblock << std::endl
				<< "cblock: " << cblock << std::endl
				<< "firstname: " << firstname << std::endl;*/
			

		}

		// Destructor
		~datasetreader() {H5Sclose(virtual_filespace);}

		// Members
		inline int rows() { return dims[0]; }
		inline int cols() { return dims[1]; }
		string filename(int i) 
		{
			//std::cout << "hello" << std::endl;
			ostringstream  file;
			file << rootname << std::setw(4) << setfill('0') << i+0 << ".h5";
			//std::cout << (nfiles > 1 ? file.str() : file.str()) << std::endl;
			return (nfiles > 1 ? file.str() : file.str());
		}

		void write(SharedMatrix<MatrixXd> M, string oname);
		void getextents(string dataset_name);
		SharedMatrix<MatrixXd> read(string dataset_name);
		void select_all();
		void select_every(int nskip);
		void select_procpart(SharedMatrix<MatrixXd> &M, hid_t *space);
	};

	void datasetreader::getextents(string dataset_name)
	{
		hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(plist, BLACS::COMM_ACTIVE, MPI::INFO_NULL);
		hid_t file = H5Fopen(firstname.c_str(), H5F_ACC_RDONLY, plist);

		hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);
		hid_t filespace = H5Dget_space(dataset);

		assert((H5Sget_simple_extent_dims(filespace, dims, NULL)==2) && "THE RANK OF THE DATASET IS NOT 2");

		H5Sclose(filespace);
		H5Dclose(dataset);
		H5Pclose(plist);
		H5Fclose(file);

		virtual_dims[0] = nfiles * dims[0];
		virtual_dims[1] = dims[1];
		virtual_filespace = H5Screate_simple(/*rank*/ 2, virtual_dims, NULL);

		
		// the snapshot matrix is savec transposed on disk, so we inverse here
		datamat.resize(virtual_dims[1], virtual_dims[0], rblock, cblock);

		// The internal matrix is transposed from row/column-major translation
		//datamat.local_matrix.resize(datamat.local_matrix.cols(), datamat.local_matrix.rows());		
	}

	void datasetreader::select_all()
	{
		

		hyperslab2D sel;		
		// Everything is transposed to match the layout on disk

		// First select the main hyperslab (full blocks)
		sel.block[0] = cblock;
		sel.block[1] = rblock;

		sel.offset[0] = cblock * BLACS::mycol;
		sel.offset[1] = rblock * BLACS::myrow;

		sel.count[0] = floor(datamat.local_matrix.cols() / cblock);
		sel.count[1] = floor(datamat.local_matrix.rows() / rblock);

		sel.stride[0] = cblock * BLACS::grid_cols;
		sel.stride[1] = rblock * BLACS::grid_rows;		

		// Only add the selection if it is non-empty
		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);

		// Second, add the remaining columns, with full row/col-block
		sel.block[0] = cblock;
		sel.block[1] = datamat.local_matrix.rows() % rblock;

		sel.offset[0] = cblock * BLACS::mycol;
		sel.offset[1] = rblock * floor(datamat.rows() / rblock);

		sel.count[0] = floor(datamat.local_matrix.cols() / cblock);
		sel.count[1] = 1;

		sel.stride[0] = cblock * BLACS::grid_cols;
		sel.stride[1] = 1;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Third, add the remaining rows, with full col/row-block
		sel.block[0] = datamat.local_matrix.cols() % cblock;
		sel.block[1] = rblock;

		sel.offset[0] = cblock * floor(datamat.cols() / cblock);
		sel.offset[1] = rblock * BLACS::myrow;

		sel.count[0] = 1;
		sel.count[1] = floor(datamat.local_matrix.rows() / rblock);

		sel.stride[0] = 1;
		sel.stride[1] = rblock * BLACS::grid_rows;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Finally, add the bottom right corner
		sel.block[0] = datamat.local_matrix.cols() % cblock;
		sel.block[1] = datamat.local_matrix.rows() % rblock;

		sel.offset[0] = cblock * floor(datamat.cols() / cblock);
		sel.offset[1] = rblock * floor(datamat.rows() / rblock);

		sel.count[0] = 1;
		sel.count[1] = 1;

		sel.stride[0] = 1;
		sel.stride[1] = 1;

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
		selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		/*hsize_t start[2], end[2];
		H5Sget_select_bounds(virtual_filespace, start, end );
		std::cout << "(" << BLACS::myrank << ")\tvirtual extents: start is " << start[0] << " x " << start[1] << ", end is " << end[0] << " x " << end[1]  << std::endl;*/
	}



	void datasetreader::select_every(int nskip)
	{

		int ccblock = nskip * cblock;


		hyperslab2D sel;		
		// Everything is transposed to match the layout on disk

		// First select the main hyperslab (full blocks)
		sel.block[0] = ccblock;
		sel.block[1] = rblock;

		sel.offset[0] = ccblock * BLACS::mycol;
		sel.offset[1] = rblock * BLACS::myrow;

		sel.count[0] = floor(datamat.local_matrix.cols() / ccblock);
		sel.count[1] = floor(datamat.local_matrix.rows() / rblock);

		sel.stride[0] = ccblock * BLACS::grid_cols;
		sel.stride[1] = rblock * BLACS::grid_rows;		

		// Only add the selection if it is non-empty
		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);

		// Second, add the remaining columns, with full row/col-block
		sel.block[0] = ccblock;
		sel.block[1] = datamat.local_matrix.rows() % rblock;

		sel.offset[0] = ccblock * BLACS::mycol;
		sel.offset[1] = rblock * floor(datamat.rows() / rblock);

		sel.count[0] = floor(datamat.local_matrix.cols() / ccblock);
		sel.count[1] = 1;

		sel.stride[0] = ccblock * BLACS::grid_cols;
		sel.stride[1] = 1;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Third, add the remaining rows, with full col/row-block
		sel.block[0] = datamat.local_matrix.cols() % ccblock;
		sel.block[1] = rblock;

		sel.offset[0] = ccblock * floor(datamat.cols() / ccblock);
		sel.offset[1] = rblock * BLACS::myrow;

		sel.count[0] = 1;
		sel.count[1] = floor(datamat.local_matrix.rows() / rblock);

		sel.stride[0] = 1;
		sel.stride[1] = rblock * BLACS::grid_rows;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Finally, add the bottom right corner
		sel.block[0] = datamat.local_matrix.cols() % ccblock;
		sel.block[1] = datamat.local_matrix.rows() % rblock;

		sel.offset[0] = ccblock * floor(datamat.cols() / ccblock);
		sel.offset[1] = rblock * floor(datamat.rows() / rblock);

		sel.count[0] = 1;
		sel.count[1] = 1;

		sel.stride[0] = 1;
		sel.stride[1] = 1;

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
		selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);


		hsize_t deb[2];
		hsize_t fin[2];
		H5Sget_select_bounds(virtual_filespace, deb, fin);

		if (BLACS::ROOT)
		{
			cout << "selection is " << deb[0] << ", " << deb[1] << " to " << fin[0] << ", " << fin[1] << endl;
		}





		// First select the individual stripes, regardless of process id
		sel.block[0] = 1;
		sel.block[1] = datamat.rows();

		sel.offset[0] = 0;
		sel.offset[1] = 0;

		sel.count[0] = ceil(datamat.cols() / nskip);
		sel.count[1] = 1;

		sel.stride[0] = nskip;
		sel.stride[1] = 1;		

		// Only add the selection if it is non-empty
		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(virtual_filespace, H5S_SELECT_AND, sel.offset, sel.stride, sel.count, sel.block);

		H5Sget_select_bounds(virtual_filespace, deb, fin);

		if (BLACS::ROOT)
		{
			cout << "selection is " << deb[0] << ", " << deb[1] << " to " << fin[0] << ", " << fin[1] << endl;
		}

		// the snapshot matrix is savec transposed on disk, so we inverse here
		datamat.resize(virtual_dims[1], ceil(virtual_dims[0] / nskip), rblock, cblock);

	}



	SharedMatrix<MatrixXd> datasetreader::read(string dataset_name)
	{
		int ncolread = 0;
		int nvalues_selected;
		hyperslab2D sel;

		sel.count[0] = datamat.local_matrix.rows() * datamat.local_matrix.cols();
		hid_t memspace = H5Screate_simple(/*rank*/ 1, sel.count, NULL);

		for (int filenum = 1; filenum <= nfiles; filenum++)
		{
			// Current file selection
			hyperslab2D sel;
			sel.block[0] = dims[0];
			sel.block[1] = dims[1];

			sel.offset[0] = (filenum-1) * dims[0];
			sel.offset[1] = 0;

			sel.count[0] = 1;
			sel.count[1] = 1;

			sel.stride[0] = 1;
			sel.stride[1] = 1;

			hid_t current_filespace = H5Scopy(virtual_filespace);
			H5Sselect_hyperslab(current_filespace, H5S_SELECT_AND, sel.offset, sel.stride, sel.count, sel.block);

			sel.offset[0] = - (filenum-1) * dims[0];
			sel.offset[1] = 0;
			if (H5Soffset_simple(current_filespace, (hssize_t*)sel.offset) < 0)
				cout << "An error happened during shift of selection" << endl << flush;

			// Current memory selection
			nvalues_selected = H5Sget_select_npoints(current_filespace);
			if (nvalues_selected < 0)
				cout << "get select npoints is negative" << endl << flush;

			sel.block[0] = nvalues_selected;

			sel.offset[0] = ncolread * datamat.local_matrix.rows();

			sel.count[0] = 1;

			sel.stride[0] = 1;

			H5Sselect_hyperslab(memspace, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);
			ncolread += nvalues_selected / datamat.local_matrix.rows();
			
			

			// Create a communicator with only processes that must read data from current file
			int iRead = (nvalues_selected != 0);
			MPI::Intracomm ReadingComm = BLACS::COMM_ACTIVE.Split(iRead, BLACS::myrank);

			
			// Only processes that must read are involved
			if (iRead)
			{
				string fname = filename(filenum);
				//cout << fname << endl << flush;
				
				hid_t plist_fa = H5Pcreate(H5P_FILE_ACCESS);
				H5Pset_fapl_mpio(plist_fa, ReadingComm, MPI::INFO_NULL);
				hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_fa);

				

				hid_t dataset = H5Dopen2(file, dataset_name.c_str(), H5P_DEFAULT);

				hid_t plist_tr = H5Pcreate(H5P_DATASET_XFER);
				H5Pset_dxpl_mpio(plist_tr, H5FD_MPIO_COLLECTIVE);

				herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, current_filespace, plist_tr, datamat.localData());

				H5Pclose(plist_tr);
				
				H5Sclose(current_filespace);
				H5Dclose(dataset);
				H5Pclose(plist_fa);
				H5Fclose(file);


				std::cout << "(" << BLACS::myrank << "|" << filenum <<")\t"<< flush;				
			}

			ReadingComm.Free();
		}

		H5Sclose(memspace);
		//datamat.local_matrix.transposeInPlace();

		return datamat;
	}

	void datasetreader::select_procpart(SharedMatrix<MatrixXd> &M, hid_t *space)
	{
		hyperslab2D sel;		

		// First select the main hyperslab (full blocks)
		sel.block[0] = M.rblock();
		sel.block[1] = M.cblock();

		sel.offset[0] = M.rblock() * BLACS::myrow;
		sel.offset[1] = M.cblock() * BLACS::mycol;

		sel.count[0] = floor(M.local_matrix.cols() / M.rblock());
		sel.count[1] = floor(M.local_matrix.rows() / M.cblock());

		sel.stride[0] = M.rblock() * BLACS::grid_rows;
		sel.stride[1] = M.cblock() * BLACS::grid_cols;		

		// Only add the selection if it is non-empty
		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(*space, H5S_SELECT_SET, sel.offset, sel.stride, sel.count, sel.block);

		// Second, add the remaining columns, with full row-block
		sel.block[0] = M.rblock();
		sel.block[1] = M.local_matrix.rows() % M.cblock();

		sel.offset[0] = M.rblock() * BLACS::myrow;
		sel.offset[1] = M.cblock() * floor(M.cols() / M.cblock());

		sel.count[0] = floor(M.local_matrix.cols() / M.rblock());
		sel.count[1] = 1;

		sel.stride[0] = M.rblock() * BLACS::grid_rows;
		sel.stride[1] = 1;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(*space, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Third, add the remaining rows, with full col-block
		sel.block[0] = M.local_matrix.cols() % M.rblock();
		sel.block[1] = M.cblock();

		sel.offset[0] = M.rblock() * floor(M.rows() / M.rblock());
		sel.offset[1] = M.cblock() * BLACS::mycol;

		sel.count[0] = 1;
		sel.count[1] = floor(M.local_matrix.rows() / M.cblock());

		sel.stride[0] = 1;
		sel.stride[1] = M.cblock() * BLACS::grid_cols;		

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
			selection.push_back(sel);*/
		H5Sselect_hyperslab(*space, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);

		// Finally, add the bottom right corner
		sel.block[0] = M.local_matrix.cols() % M.rblock();
		sel.block[1] = M.local_matrix.rows() % M.cblock();

		sel.offset[0] = M.rblock() * floor(M.rows() / M.rblock());
		sel.offset[1] = M.cblock() * floor(M.cols() / M.cblock());

		sel.count[0] = 1;
		sel.count[1] = 1;

		sel.stride[0] = 1;
		sel.stride[1] = 1;

		/*if (sel.block[0]*sel.block[1] * sel.count[0]*sel.count[1])
		selection.push_back(sel);*/
		H5Sselect_hyperslab(*space, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);
	}

	void datasetreader::write(SharedMatrix<MatrixXd> M, string oname)
	{
		hyperslab2D sel;
		
		sel.count[0] = M.rows();
		sel.count[1] = M.cols();
		hid_t filespace = H5Screate_simple(/*rank*/ 2, sel.count, NULL);
		
		int iWrite = (H5Sget_select_npoints(filespace) != 0);
		MPI::Intracomm WritingComm = BLACS::COMM_ACTIVE.Split(iWrite, BLACS::myrank);

		if (iWrite)
		{
			std::cout << "(" << BLACS::myrank << ")\t iWrite"<< flush;
			M.local_matrix.transposeInPlace();
			std::cout << "(" << BLACS::myrank << ")\ttransposing done"<< flush;

			hid_t plist_fa = H5Pcreate(H5P_FILE_ACCESS);
			H5Pset_fapl_mpio(plist_fa, WritingComm, MPI::INFO_NULL);
			hid_t file = H5Fcreate(oname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_fa);

			hid_t dataset = H5Dcreate(file, "snapshots", H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
						
			sel.count[0] = M.local_matrix.rows();
			sel.count[1] = M.local_matrix.cols();
			hid_t memspace = H5Screate_simple(/*rank*/ 2, sel.count, NULL);
			H5Sselect_all(memspace);

			select_procpart(M, &filespace);

			hid_t plist_tr = H5Pcreate(H5P_DATASET_XFER);
			H5Pset_dxpl_mpio(plist_tr, H5FD_MPIO_COLLECTIVE);

			std::cout << "(" << BLACS::myrank << ")\twriting"<< flush;
			H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_tr, M.localData());
			std::cout << "(" << BLACS::myrank << ")\twriting done"<< flush;

			H5Pclose(plist_tr);
			H5Sclose(memspace);
			H5Dclose(dataset);
			H5Fclose(file);
			H5Pclose(plist_fa);
		}
		H5Sclose(filespace);
		WritingComm.Free();
	}

}	// end namespace phdfp
#endif // PHDFP_DATASET