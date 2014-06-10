// DMD.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"



#ifndef USE_PRECOMPILED_HEADER
#include "mpi.h"

#include "boost/lexical_cast.hpp"

#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include "PBLAS.h"
#include "BLACS.h"

#include "H5inDMD.h"

#define USE_PROFILER
#define USE_BOOST_CHRONO
#include "tic-toc-profiler.hpp"

#endif



#include "SharedMatrix.h"
#include "ScaSVD.h"

#include "ScaSolve.h"

#include "ScaEigenSolver.h"



using namespace std;
using namespace Eigen;
using namespace peigen;
using namespace phdfp;


int main(int argc, char* argv[])
{
	int mkl_res = mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
	MPI::Init();
	int rank, numtasks;
	rank = MPI::COMM_WORLD.Get_rank();
	numtasks = MPI::COMM_WORLD.Get_size();
	bool ROOT = (rank == 0);

	DoProfiler prof;

	cout.precision(2*std::numeric_limits< double >::digits10);

	// Create the BLACS grid
	BLACS::init(numtasks);
	if (BLACS::ROOT)
		std::cout << "BLACS Initialized. MKL says " << mkl_res << endl << flush;
	//BLACS::printGrid();

	MPI::COMM_WORLD.Barrier(); // For printing purposes


	// Input arguments
	const int nfiles = 7;// boost::lexical_cast<int>(argv[1]);
	const int nskip_step = 3;// boost::lexical_cast<int>(argv[2]);

	if (BLACS::ROOT)
	{
		std::cout << "Input arguments: " << nfiles << ", " << nskip_step << endl << flush;
	}



	// Use all the processes for faster IO, regardless of them being used in the process grid
	prof.tic("Read");
	if (BLACS::ROOT)
		std::cout << "Creating Reader" << endl << flush;

	datasetreader dreader(nfiles, "D:/DMD/testmat");
	if (BLACS::ROOT)
		std::cout << "Datareader created." << endl << flush;

	MPI::COMM_WORLD.Barrier(); // For printing purposes

	dreader.read("testgroup_T");

	MPI::COMM_WORLD.Barrier(); // For printing purposes

	SharedMatrix<MatrixXd> snaps(dreader.createShared(6, 6, nskip_step));

	cout << snaps << endl;

	//std::cout << "From the return value: " << prof.toc("Read") << endl << flush;
	prof.toc("Read", "The time Reading took is: ");



	// Discard inactive processes immediately, this is to avoid crashes caused e.g. by inactive process calling barrier(ctxt, All), or descinit()
	// Note: If needed, use the MPI::Intracomm BLACS::COMM_ACTIVE to avoid deadlocks with incative processes
	if (BLACS::active) // Only processes that have a place in the grid
	{

		//

		//dreader.getextents("snapshots_T");

		////cout << "rows: " << dreader.rows() << ", cols: " << dreader.cols() << endl << flush;

		//dreader.select_every(nskip_step);
		////dreader.select_all();
		//dreader.read("snapshots_T");

		////std::cout << "READING DONE." << endl << flush;


		//SharedMatrix<MatrixXd> snaps = dreader.datamat;
		//prof.toc("Read");

		cout << snaps;

		///**************************************************************************************************/
		///*-------------------------------------     START OF DMD     -------------------------------------*/
		///**************************************************************************************************/


		int Np = snaps.rows() / 4;
		int Nt = snaps.cols();
		BLACS::COMM_ACTIVE.Barrier();

		/////**************************************************************************************************/
		/////*-------------------------------------       DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/
		prof.tic("SVD");
		ScaSVD<MatrixXd> svd(snaps.block(0, 0, 4 * Np, Nt - 1), true, true);
		snaps.clear();

		if (BLACS::myrank == 0)
			cout << "SVD done" << endl << flush;


		BLACS::COMM_ACTIVE.Barrier();
		if (ROOT)
			cout << "singular values on process  " << rank << " are " << svd.singularValues.transpose() << endl << flush;


		if (ROOT)
			cout << "HERE COMES U" << endl << flush;
		cout << svd.matrixU;
		

		if (ROOT)
			cout << "HERE COMES Vt" << endl << flush;
		cout << svd.matrixVt;


		double r_svd = svd.residual(snaps.block(0, 0, 4 * Np, Nt - 1));
		//if (ROOT)
		cout << "Residual from SVD: " << r_svd << endl << flush;

		

		prof.toc("SVD");
		/////**************************************************************************************************/
		/////*-------------------------------------      /DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/

		BLACS::COMM_ACTIVE.Barrier();


		///////**************************************************************************************************/
		///////*-----------------------------      DO B := Ut * S2 * V * SIG+     ------------------------------*/
		///////**************************************************************************************************/
		prof.tic("MultiplyB");
		//SharedMatrix<MatrixXd> S2(4*Np, NtperF - 1);
		// Do a virtual shift
		/*S2.local_matrix = MatrixXd::Ones(S2.local_matrix.rows(), S2.local_matrix.cols()) * rank;
		for (int c = S2.cblock-1; c < S2.local_matrix.cols(); c+=S2.cblock)
		{
		S2.local_matrix.col(c) = MatrixXd::Ones(S2.local_matrix.rows(), 1) * (BLACS::grid_cols*BLACS::myrow+((rank+1)%BLACS::grid_cols));
		}
		if (ROOT)


		cout << "HERE COMES S2" << endl << flush;
		cout << S2 ;*/


		//cout << "HERE COMES U BEFORE" << endl << flush;
		//cout << svd.matrixU << endl << endl;
		SharedMatrix<MatrixXd> tmpMat = svd.matrixU.transpose() * snaps.block(0, 1, 4 * Np, Nt - 1);
		
		cout << "HERE COMES tmpMat" << endl << flush;
		cout << tmpMat << endl << endl;
		
		
		snaps.clear();
		SharedMatrix<MatrixXd> B = tmpMat * svd.matrixVt.transpose();



		svd.matrixU.clear();

		if (ROOT)
			cout << "HERE COMES Partial B " << endl << flush;
		cout << B;


		//// Now only a right-multiply by SIG+ is left to do
		//// Each process owns the entire vector sig. Figuring out which value goes to which process
		//// and creating a shared matrix and multiplying is a pain. 
		//// Since SIG+ is diagonal we do what it does: scale the columns of B by the propre value.
		MatrixXd SIGplus = MatrixXd::Zero(B.local_matrix.cols(), B.local_matrix.cols());
		double pinv_tol = std::numeric_limits<double>::epsilon() * (4 * Np) * svd.singularValues(0);
		if (ROOT)
		{
			cout << " pinv_tol: " << pinv_tol << endl << flush;
			//assert(pinv_tol == 1.4825790295167806e-011 && "Tolerance for pseudo-inverse differs");
		}
		for (int i = 0; i < B.local_matrix.cols(); i++)
		{
			int index = i % B.cblock() + (floor(i / B.cblock())*BLACS::grid_cols + BLACS::mycol) * B.cblock();
			if (svd.singularValues(index) > pinv_tol)
				SIGplus(i, i) = 1 / svd.singularValues(index);
		}

		if (ROOT)
		{
			cout << "HERE COMES SIG+" << endl << flush;
			cout << SIGplus.diagonal().transpose() << endl;
		}


		B.local_matrix = B.local_matrix * SIGplus;
		if (ROOT)
			cout << "B is computed " << endl;
		if (ROOT)
			cout << "HERE COMES Ut M V SIG+" << endl << flush;
		cout << B;

		prof.toc("MultiplyB");
		/////**************************************************************************************************/
		/////*------------------------------     /DO B := Ut * M * V * SIG+     ------------------------------*/
		/////**************************************************************************************************/


		BLACS::COMM_ACTIVE.Barrier();

		/////**************************************************************************************************/
		/////*------------------------------     DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/
		prof.tic("EigenProblem");
		//ScaEigenSolve<MatrixXd> pes(B);

		

		ScaEigenSolver<MatrixXd> eig(B, true, EigSchur);

		double r_eig = eig.global_residual(B);
		if (ROOT)
			cout << "Residual from Eigen problem: " << r_eig << endl << flush;

		// Matrix of eigen vectors
		SharedMatrix<MatrixXcd> X = eig.eigenVectors();
		MatrixXcd lambdas = eig.eigenValues(); 

		prof.toc("EigenProblem");
		/////**************************************************************************************************/
		/////*------------------------------    /DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/

		BLACS::COMM_ACTIVE.Barrier();


		/////**************************************************************************************************/
		/////*------------------------------      DO A LINEAR SYSTEM SOLVE      ------------------------------*/
		/////**************************************************************************************************/
		prof.tic("LinearSolve");
		//cout << "(" << BLACS::myrank << ")" << endl;
		SharedMatrix<MatrixXd> rhs = svd.matrixU.transpose() * snaps.block(0, Nt - 1, 4 * Np, 1);
		svd.matrixU.clear();
		snaps.clear();
		SharedMatrix<MatrixXcd> rhsZ = rhs.cast<std::complex<double> >();

		BLACS::COMM_ACTIVE.Barrier();

		ScaSolve<MatrixXcd> solver(X, rhsZ, peigen::pxgesv);

		if (ROOT)
		{
			cout << "HERE COME the solution" << endl << flush;
		}
		cout << solver.solution;
		
		cout << "Residual from the system solve: " << solver.residual(X, rhsZ) << endl << flush;

		/////**************************************************************************************************/
		/////*------------------------------      /DO A LINEAR SYSTEM SOLVE      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*-----------------------------      APPLY WEIGHT TO THE MODES      ------------------------------*/
		/////**************************************************************************************************/

		SharedMatrix<MatrixXcd> Modes = svd.matrixU.cast<std::complex<double> >() * X;

		Modes.ColScale(solver.solution);

		
		//cout << Modes;
		prof.toc("LinearSolve");
		/////**************************************************************************************************/
		/////*-----------------------------      /APPLY WEIGHT TO THE MODES      -----------------------------*/
		/////**************************************************************************************************/

		//BLACS::COMM_ACTIVE.Barrier();


		/////**************************************************************************************************/
		/////*------------------------------      Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/
		//prof.tic("Energy");
		//SharedMatrix<MatrixXd> amplitudes(1, Modes.cols(), 1, Modes.cblock());

		//amplitudes.local_matrix = Modes.local_matrix.colwise().squaredNorm();

		//char scope[7] = { 'C', 'O', 'L', 'U', 'M', 'N', '\0' };
		//char top[2] = { ' ', '\0' };
		//// Sum the norm of columns by process column
		//BLACS::Cdgsum2d(BLACS::ctxt, scope, top, 1, amplitudes.local_matrix.cols(), amplitudes.localData(), 1, -1, -1);

		//// cout the amplitudes
		//for (int rr = 0; rr<BLACS::grid_rows; rr++)
		//{
		//	for (int cc = 0; cc<BLACS::grid_cols; cc++)
		//	{
		//		//if((BLACS::myrow==rr) && (BLACS::mycol==cc))
		//		//cout << "(" << BLACS::myrow << ","<<BLACS::mycol << ") " << amplitudes.local_matrix << flush;
		//		BLACS::COMM_ACTIVE.Barrier();
		//	}
		//	//if(BLACS::myrank==0)
		//	//cout << endl << flush;
		//	BLACS::COMM_ACTIVE.Barrier();
		//}

		//// Gather the global amplitudes vector, on a row basis
		//amplitudes.global_matrix.resize(1, amplitudes.cols());
		//int ncols = Modes.cols();
		//int modes_cblock = Modes.cblock();
		//int offsetmpi = 0;
		//for (int cc = 0; cc<BLACS::grid_cols; cc++)
		//{
		//	int n = BLACS::numroc_(&ncols, &modes_cblock, &cc, BLACS::iZERO, &(BLACS::grid_cols));
		//	int pr = BLACS::Cblacs_pnum(BLACS::ctxt, BLACS::myrow, cc);
		//	BLACS::COMM_ACTIVE.Irecv(amplitudes.globalData() + offsetmpi, n, MPI::DOUBLE, pr, pr);
		//	//cout << BLACS::myrank <<" Irecv from " << pr << " for " << n << " values with offset"<<offsetmpi<<endl << flush;
		//	offsetmpi += n;
		//}


		//for (int cc = 0; cc<BLACS::grid_cols; cc++)
		//{
		//	int pr = BLACS::Cblacs_pnum(BLACS::ctxt, BLACS::myrow, cc);
		//	BLACS::COMM_ACTIVE.Send(amplitudes.localData(), amplitudes.local_matrix.cols(), MPI::DOUBLE, pr, BLACS::myrank);
		//	//cout << BLACS::myrank <<" send " << ampl << " to " << pr << " for " << ampl.cols() << " values"<<endl << flush;
		//}

		///*for (int proc = 0; proc < BLACS::grid_rows*BLACS::grid_cols; proc++)
		//{
		//amplitudes.gather(proc);
		//}*/

		////BLACS::COMM_ACTIVE.Barrier();

		//if (ROOT)
		//	cout << "amplitudes are done" << endl << flush;

		///*for(int rr = 0; rr<BLACS::grid_rows*BLACS::grid_cols; rr++)*/
		//for (int rr = 0; rr<1; rr++)
		//{
		//	if (BLACS::myrank == rr)
		//		cout << "(" << BLACS::myrank << ") " << amplitudes.global_matrix.transpose() << endl << flush;
		//	BLACS::COMM_ACTIVE.Barrier();
		//}
		//prof.toc("Energy");
		/////**************************************************************************************************/
		/////*-----------------------------      /Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/



		//BLACS::COMM_ACTIVE.Barrier();



		/////**************************************************************************************************/
		/////*------------------------------       PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/
		//prof.tic("DumpH5");
		//string out_file = "DMDosc.h5";
		//BLACS::COMM_ACTIVE.Barrier();
		//// Find which column of Modes has the most energy
		//int NMODES = 3;
		//MatrixXd::Index i_mode;

		//if (ROOT)
		//{
		//	cout << "(" << BLACS::myrank << ")" << endl;
		//	Matrix<double, Dynamic, 2, RowMajor> spectrum(Nt - 1, 2);
		//	spectrum.col(0) = lambdas.real().cwiseQuotient(lambdas.cwiseAbs()).array().acos().matrix();
		//	spectrum.col(1) = amplitudes.global_matrix.row(0);

		//	std::ofstream s("spectrum.txt");
		//	if (s.is_open())
		//	{
		//		s << spectrum << '\n';
		//		s.close();
		//	}

		//	std::ofstream l("lambdas.txt");
		//	if (l.is_open())
		//	{
		//		l << lambdas << '\n';
		//		l.close();
		//	}

		//	hyperslab2D sel;
		//	sel.count[0] = spectrum.rows();
		//	sel.count[1] = spectrum.cols();
		//	//hid_t filespace = H5Screate_simple(2, sel.count, NULL); 
		//	//hid_t memspace = H5Screate_simple(/*rank*/ 2, sel.count, NULL);

		//	hid_t file = H5Fcreate(out_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		//	//	hid_t dataset = H5Dcreate(file, "/spectrum", H5T_NATIVE_DOUBLE, filespace,
		//	//				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//	//	herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
		//	//				H5P_DEFAULT, spectrum.data());

		//	//	H5Dclose(dataset);
		//	H5Fclose(file);/*
		//				   H5Sclose(memspace);
		//				   H5Sclose(filespace);*/
		//}

		//for (int mode = 0; mode < NMODES; mode++)
		//{
		//	amplitudes.global_matrix.row(0).maxCoeff(&i_mode);

		//	while (lambdas(i_mode).imag() < 0)
		//	{
		//		amplitudes.global_matrix(i_mode) = -1;
		//		amplitudes.global_matrix.row(0).maxCoeff(&i_mode);
		//	}
		//	amplitudes.global_matrix(i_mode) = -1;

		//	if (lambdas(i_mode).imag() >= 0)
		//	{
		//		//cout << rank<< " i_mode " << i_mode << endl <<flush;

		//		MPI::Intracomm ColComm = BLACS::COMM_ACTIVE.Split(BLACS::mycol, BLACS::myrow);
		//		int icol = (i_mode / Modes.cblock()) % BLACS::grid_cols;
		//		int i_loc = (i_mode / (Modes.cblock() * BLACS::grid_cols)) + (i_mode % Modes.cblock());
		//		//cout << rank<< " icol " << icol << endl <<flush;
		//		//cout << rank<< " i_loc " << i_loc << endl <<flush;

		//		if (BLACS::mycol == icol)
		//		{
		//			hyperslab2D sel;
		//			//cout << rank << " in "  << endl <<flush;
		//			/*
		//			* Set up file access property list with parallel I/O access
		//			*/
		//			hid_t plist_fa = H5Pcreate(H5P_FILE_ACCESS);
		//			H5Pset_fapl_mpio(plist_fa, ColComm, MPI::INFO_NULL);
		//			//cout << rank<< " before create "  << endl <<flush;
		//			/*
		//			* Create a new file collectively and release property list identifier.
		//			*/
		//			hid_t file = H5Fopen(out_file.c_str(), H5F_ACC_RDWR, plist_fa);
		//			H5Pclose(plist_fa);
		//			//cout << rank<< " after create "  << endl <<flush;
		//			/*
		//			* Create the dataspace for the dataset.
		//			*/
		//			sel.count[0] = 4 * Np;
		//			sel.count[1] = 2;
		//			hid_t filespace = H5Screate_simple(2, sel.count, NULL);

		//			/*
		//			* Create the dataset with default properties and close filespace.
		//			*/
		//			ostringstream  datasetname;
		//			datasetname << "/mode" << mode;

		//			hid_t dataset = H5Dcreate(file, datasetname.str().c_str(), H5T_NATIVE_DOUBLE, filespace,
		//				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//			//H5Sclose(filespace);
		//			//cout << rank<< " dataset created "  << endl <<flush;
		//			/*
		//			* Each process defines dataset in memory and writes it to the hyperslab
		//			* in the file.
		//			*/
		//			sel.count[0] = Modes.local_matrix.rows();
		//			sel.count[1] = 2;
		//			sel.offset[0] = 0;
		//			sel.offset[1] = 0;
		//			hid_t memspace = H5Screate_simple(/*rank*/ 2, sel.count, NULL);
		//			H5Sselect_hyperslab(memspace, H5S_SELECT_SET, sel.offset, NULL, sel.count, NULL);
		//			ColComm.Barrier();
		//			/*H5Sget_select_bounds(memspace, start, end );
		//			cout << "loop " << MPI::COMM_WORLD.Get_rank() << ": start is " << start[0] << "x" << start[1] << ", end is " << end[0] << "x" << end[1] <<endl;
		//			ColComm.Barrier();*/

		//			/*
		//			* Select hyperslab in the file.
		//			*/

		//			sel.count[0] = 1;
		//			sel.count[1] = 2;
		//			sel.offset[0] = Modes.rblock() * (BLACS::myrow);
		//			sel.offset[1] = 0;
		//			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, sel.offset, NULL, sel.count, NULL);


		//			for (int r = 0; r < min(Modes.rblock(), (int)Modes.local_matrix.rows()); r++)
		//			{
		//				sel.stride[0] = Modes.rblock() * BLACS::grid_rows;
		//				sel.stride[1] = 1;
		//				sel.count[0] = ceil((double)(Modes.local_matrix.rows() - r) / Modes.rblock());
		//				sel.count[1] = 2;
		//				sel.offset[0] = r + Modes.rblock() * BLACS::myrow;
		//				sel.offset[1] = 0;
		//				sel.block[0] = 1;
		//				sel.block[1] = 1;
		//				H5Sselect_hyperslab(filespace, H5S_SELECT_OR, sel.offset, sel.stride, sel.count, sel.block);
		//			}
		//			ColComm.Barrier();
		//			/*H5Sget_select_bounds(filespace, start, end );
		//			cout << "filespace " << MPI::COMM_WORLD.Get_rank() << ": start is " << start[0] << "x" << start[1] << ", end is " << end[0] << "x" << end[1] <<endl;
		//			ColComm.Barrier();*/

		//			/*
		//			* Create property list for collective dataset write.
		//			*/
		//			MatrixXd CurrentMode(Modes.local_matrix.rows(), 2);
		//			CurrentMode.col(0) = Modes.local_matrix.col(i_loc).cwiseAbs();
		//			CurrentMode.col(1) = Modes.local_matrix.col(i_loc).imag().binaryExpr(Modes.local_matrix.col(i_loc).real(), std::ptr_fun(atan2<double, double>));

		//			/*if (ROOT)
		//			cout << CurrentMode.block(0,0, 40, 2);*/
		//			if (ROOT)
		//				cout << "Writing mode " << mode << "..." << endl;

		//			CurrentMode.transposeInPlace();



		//			hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
		//			H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		//			//cout << "Just before write " <<endl << flush;
		//			herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
		//				plist_id, CurrentMode.data());
		//			//cout << "Just after write " <<endl << flush;
		//			/*
		//			* Close/release resources.
		//			*/
		//			H5Pclose(plist_id);
		//			H5Sclose(memspace);
		//			H5Dclose(dataset);
		//			H5Sclose(filespace);

		//			//H5Pclose(plist_id);
		//			H5Fclose(file);
		//		}
		//	}
		//}
		//prof.toc("DumpH5");
		/////**************************************************************************************************/
		/////*------------------------------      /PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/

		stringstream profile_data;
		profile_data << argv[0] << "-" << rank << ".yml";

		//if (BLACS::myrank == 0)
		prof.dump(profile_data.str());

		BLACS::COMM_ACTIVE.Barrier();
	} // end IF ACTIVE
	else
	{
		cout << "Pfff" << endl;
	}

	if (ROOT)	// not having if(root) prevents segfault at the begining of the program (!?!)
		cout << "DONE" << endl;

	BLACS::finalize();
	MPI::Finalize();

	return 0;
}

