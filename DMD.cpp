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

#include "ColumnSquaredNorm.h"
#endif



#include "SharedMatrix.h"
#include "ScaSVD.h"

#include "ScaSolve.h"

#include "ScaEigenSolver.h"
#include "Vandermonde.h"


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
	const int nfiles = 1;// boost::lexical_cast<int>(argv[1]);
	const int nskip_step = 5;// boost::lexical_cast<int>(argv[2]);

	if (BLACS::ROOT)
	{
		std::cout << "Input arguments: " << nfiles << ", " << nskip_step << endl << flush;
	}



	// Use all the processes for faster IO, regardless of them being used in the process grid
	prof.tic("Read");
	if (BLACS::ROOT)
		std::cout << "Creating Reader" << endl << flush;

	datasetreader dreader(nfiles, "D:/DMD/DMD/x64/NNDEB/Re350_oscillating");
	if (BLACS::ROOT)
		std::cout << "Datareader created." << endl << flush;

	MPI::COMM_WORLD.Barrier(); // For printing purposes

	dreader.read("snapshots_T");

	MPI::COMM_WORLD.Barrier(); // For printing purposes

	SharedMatrix<MatrixXd> snaps(dreader.createShared(6, 6, nskip_step));

	//cout << snaps << endl;

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

		//cout << snaps;

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


		/*if (ROOT)
			cout << "HERE COMES U" << endl << flush;
		cout << svd.matrixU;*/
		

		/*if (ROOT)
			cout << "HERE COMES Vt" << endl << flush;
		cout << svd.matrixVt;*/


		double r_svd = svd.residual(snaps.block(0, 0, 4 * Np, Nt - 1));
		//if (ROOT)
		cout << "Residual from SVD: " << r_svd << endl << flush;

		cout << "Residual GLOBAL from SVD: " << svd.global_residual(snaps.block(0, 0, 4 * Np, Nt - 1)) << endl << flush;

		
		prof.toc("SVD", "SVD done in: ");
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
		
		/*cout << "HERE COMES tmpMat" << endl << flush;
		cout << tmpMat << endl << endl;*/
		
		
		snaps.clear();
		SharedMatrix<MatrixXd> B = tmpMat * svd.matrixVt.transpose();



		svd.matrixU.clear();

		/*if (ROOT)
			cout << "HERE COMES Partial B " << endl << flush;
		cout << B;*/


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

		/*if (ROOT)
		{
			cout << "HERE COMES SIG+" << endl << flush;
			cout << SIGplus.diagonal().transpose() << endl;
		}*/


		B.local_matrix = B.local_matrix * SIGplus;
		if (ROOT)
			cout << "B is computed " << endl;
		/*if (ROOT)
			cout << "HERE COMES Ut M V SIG+" << endl << flush;
		cout << B;*/

		prof.toc("MultiplyB", "B created in: ");
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

		prof.toc("EigenProblem", "Eigen problem solved in: ");
		/////**************************************************************************************************/
		/////*------------------------------    /DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/

		BLACS::COMM_ACTIVE.Barrier();


		/////**************************************************************************************************/
		/////*------------------------------      DO A LINEAR SYSTEM SOLVE      ------------------------------*/
		/////**************************************************************************************************/
		prof.tic("LinearSolve");
		//cout << "(" << BLACS::myrank << ")" << endl;
		SharedMatrix<MatrixXd> rhs = svd.matrixU.transpose() * snaps.block(0, 0/*Nt - 1*/, 4 * Np, 1);
		svd.matrixU.clear();
		snaps.clear();
		SharedMatrix<MatrixXcd> rhsZ = rhs.cast<std::complex<double> >();

		
		
		// Construct a system so that the weights will have to be in complex conjugate pairs
		SharedMatrix<MatrixXd> System(X.rows(), X.cols(), X.rblock(), X.cblock());
		for (int k = 0; k < lambdas.rows(); ++k)
		{
			// First item of a conjugate pair
			if (lambdas(k, 0).imag() != 0)
			{
				if (BLACS::indxg2p(k, System.cblock(), BLACS::grid_cols) == BLACS::mycol)
				{
					int l = BLACS::indxg2l(k, System.cblock(), BLACS::grid_cols);
					System.local_matrix.col(l) = 2 * X.local_matrix.col(l).real();
				}
				++k;
				// Second item of a conjugate pair
				if (BLACS::indxg2p(k, System.cblock(), BLACS::grid_cols) == BLACS::mycol)
				{
					int l = BLACS::indxg2l(k, System.cblock(), BLACS::grid_cols);
					System.local_matrix.col(l) = 2 * X.local_matrix.col(l).imag();
				}
			}
			else // real eigenvalue
			{
				if (BLACS::indxg2p(k, System.cblock(), BLACS::grid_cols) == BLACS::mycol)
				{
					int l = BLACS::indxg2l(k, System.cblock(), BLACS::grid_cols);
					System.local_matrix.col(l) = X.local_matrix.col(l).real();
				}
			}
		}




		BLACS::COMM_ACTIVE.Barrier();

		//ScaSolve<MatrixXcd> solver(X, rhsZ, peigen::EigenSVD);

		ScaSolve<MatrixXd> solver(System, rhs, peigen::pxgesvx);
		cout << "Residual from the system solve: " << solver.residual(System, rhs) << endl << flush;

		// Reconstitute the solution to the original system
		SharedMatrix<MatrixXcd> weights(solver.solution.rows(), solver.solution.cols(), solver.solution.rblock(), solver.solution.cblock());
		for (int k = 0; k < lambdas.rows(); ++k)
		{
			// First item of a conjugate pair
			if (lambdas(k, 0).imag() != 0)
			{
				if (BLACS::indxg2p(k, System.rblock(), BLACS::grid_rows) == BLACS::myrow)
				{
					if (BLACS::indxg2p(k + 1, System.rblock(), BLACS::grid_rows) == BLACS::myrow)
					{
						int l = BLACS::indxg2l(k, System.rblock(), BLACS::grid_rows);
						weights.local_matrix(l, 0) = complex<double>(solver.solution.local_matrix(l, 0), solver.solution.local_matrix(l + 1, 0));
						weights.local_matrix(l + 1, 0) = complex<double>(solver.solution.local_matrix(l, 0), -solver.solution.local_matrix(l + 1, 0));
					}
					else
					{
						int ownernext = BLACS::indxg2p(k + 1, System.rblock(), BLACS::grid_rows);
						double im;
						BLACS::COMM_ACTIVE.Recv(&im, 1, MPI::DOUBLE, ownernext, ownernext);

						int l = BLACS::indxg2l(k, System.rblock(), BLACS::grid_rows);
						weights.local_matrix(l, 0) = complex<double>(solver.solution.local_matrix(l, 0), im);
					}
				}
				else
				{
					if (BLACS::indxg2p(k + 1, System.rblock(), BLACS::grid_rows) == BLACS::myrow)
					{
						int owner = BLACS::indxg2p(k + 1, System.rblock(), BLACS::grid_rows);
						double re;
						BLACS::COMM_ACTIVE.Recv(&re, 1, MPI::DOUBLE, owner, owner);

						int l = BLACS::indxg2l(k + 1, System.rblock(), BLACS::grid_rows);
						weights.local_matrix(l + 1, 0) = complex<double>(re, solver.solution.local_matrix(l + 1, 0));
					}
				}
				++k;
			}
			else // real eigenvalue
			{
				if (BLACS::indxg2p(k, System.cblock(), BLACS::grid_cols) == BLACS::mycol)
				{
					int l = BLACS::indxg2l(k, System.rblock(), BLACS::grid_rows);
					weights.local_matrix(l, 0) = complex<double>(solver.solution.local_matrix(l, 0), 0);
				}
			}
		}

		if (ROOT)
		{
			cout << "HERE COME the solution" << endl << flush;
		}
		cout << weights << endl;
		
		

		/////**************************************************************************************************/
		/////*------------------------------      /DO A LINEAR SYSTEM SOLVE      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*-----------------------------      APPLY WEIGHT TO THE MODES      ------------------------------*/
		/////**************************************************************************************************/

		SharedMatrix<MatrixXcd> Modes = svd.matrixU.cast<std::complex<double> >() * X;

		Modes.ColScale(weights);


		SharedMatrix<MatrixXcd> Vandermonde = vander<MatrixXcd>(lambdas, Modes.cols(), Modes.rblock(), Modes.cblock());
		//cout << Vandermonde << endl;

		SharedMatrix<MatrixXcd> reconstruct = snaps.cast<complex<double>>().block(0,0, 4*Np, Nt-1);
		reconstruct.pgemm(1., Modes, Vandermonde, -1.);
		//cout << reconstruct << endl;
		cout << "Residual from Modes: " << reconstruct.localBlock().cwiseAbs().maxCoeff() << endl;
		
		//cout << Modes;
		prof.toc("LinearSolve", "Eigen problem solved in: ");
		/////**************************************************************************************************/
		/////*-----------------------------      /APPLY WEIGHT TO THE MODES      -----------------------------*/
		/////**************************************************************************************************/

		//BLACS::COMM_ACTIVE.Barrier();


		/////**************************************************************************************************/
		/////*------------------------------      Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/
		prof.tic("Energy");
		MatrixXd amplitudes = ColumnSquaredNorm(Modes);
		prof.toc("Energy");

		if (BLACS::myrank == 0)
			cout << amplitudes << endl << flush;
		
		/////**************************************************************************************************/
		/////*-----------------------------      /Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*----------------------------------      Print light data      ----------------------------------*/
		/////**************************************************************************************************/
		
		if (ROOT)
		{
			cout << "Printing spectrum...\t";
			Matrix<double, Dynamic, 2, RowMajor> spectrum(Nt - 1, 2);
			spectrum.col(0) = lambdas.real().cwiseQuotient(lambdas.cwiseAbs()).array().acos().matrix();
			spectrum.col(1) = amplitudes.transpose();

			std::ofstream s("spectrum.txt");
			if (s.is_open())
			{
				s << spectrum << '\n';
				s.close();
			}
			cout << "DONE." << endl;

			cout << "Printing eigenvalues...\t";
			std::ofstream l("eigenvalues.txt");
			if (l.is_open())
			{
				l << lambdas << '\n';
				l.close();
			}
			cout << "DONE." << endl;
		}

		/////**************************************************************************************************/
		/////*----------------------------------      /Print light data      ---------------------------------*/
		/////**************************************************************************************************/
		
		//BLACS::COMM_ACTIVE.Barrier();



		/////**************************************************************************************************/
		/////*------------------------------       PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/
		prof.tic("DumpModes");
		BLACS::COMM_ACTIVE.Barrier();
		// Find which column of Modes has the most energy
		int NMODES = 1; // Number of modes to print
		MatrixXd::Index i_mode, x_mode;

		stringstream variables;

		for (int m = 0; m < NMODES; ++m)
		{
			cout << "writing mode " << m << "...";
			//Find mode with highest amplitude
			amplitudes.maxCoeff(&x_mode, &i_mode); //x_mode is always 0 since amplitudes is a 1xN matrix
			
			// Only print 1 mode out of a possible conjugate pair
			while (lambdas(i_mode, 0).imag() < 0)
			{
				amplitudes(0, i_mode) = -1;	// Get the bottom part of the spectrum out of the game
				amplitudes.maxCoeff(&x_mode, &i_mode);
			}
			amplitudes(0, i_mode) = -1;

			stringstream filename;
			filename << "mode" << m << ".abs";
			FILE *pFile;
			pFile = fopen(filename.str().c_str(), "wb");
			

			char text[81];
			sprintf(text, "Module of Mode %06i%-58s\n", m, "");
			fwrite(text, 1, 80*sizeof(char), pFile);

			sprintf(text, "%-79s\n", "part");
			fwrite(text, 1, 80 * sizeof(char), pFile);

			const int part_number = 1;
			fwrite(&part_number, 1, 1 * sizeof(int), pFile);

			sprintf(text, "%-79s\n", "hexa8");
			fwrite(text, 1, 80 * sizeof(char), pFile);

			//cout << "test 1" << endl;

			MatrixXf export = Modes.local_matrix.col(i_mode).cwiseAbs().cast<float>();
			fwrite(export.data(), 1, Modes.rows() * sizeof(float) / 4 , pFile);

			//cout << "test 2" << endl;

			fclose(pFile);
			variables << "scalar per element: Mode" << m << "abs " << filename.str() << endl;
			//cout << " abs" << m;

			/*stringstream filenameIM;
			filenameIM << "mode" << m << ".ang";
			
			pFile = fopen(filenameIM.str().c_str(), "wb");

			sprintf(text, "Angle of Mode %06i %-58s\n", m, "");
			fwrite(text, 1, 80 * sizeof(char), pFile);

			sprintf(text, "%-79s\n", "part");
			fwrite(text, 1, 80 * sizeof(char), pFile);

			fwrite(&part_number, 1, 1 * sizeof(int), pFile);

			sprintf(text, "%-79s\n", "hexa8");
			fwrite(text, 1, 80 * sizeof(char), pFile);

			export = Modes.local_matrix.col(i_mode).imag().binaryExpr(Modes.local_matrix.col(i_mode).real(), std::ptr_fun(atan2<double, double>)).cast<float>();
			fwrite(export.data(), 1, Modes.rows() * sizeof(float) / 4, pFile);

			fclose(pFile);
			variables << "scalar per element: Mode" << m << "ang " << filenameIM.str() << endl;*/

			cout << "\tDONE" << endl;
		}

		std::ofstream ofs("dmd.case", std::ofstream::out);

		ofs << "FORMAT" << endl
			<< "type: ensight gold" << endl
			<< "GEOMETRY" << endl
			<< "model: dmd.geo" << endl
			<< "VARIABLE" << endl
			<< variables.str()
			<< "TIME" << endl
			<< "time set: 1 \nnumber of steps: 1 \nfilename start number: 0 \nfilename increment: 1 \ntime values: \n0" << endl;


		ofs.close();

		prof.toc("DumpModes");
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
		cout << "PROGRAM FINISHED" << endl;

	BLACS::finalize();
	MPI::Finalize();

	return 0;
}

