// DMD.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"

#ifndef USE_PRECOMPILED_HEADER

// External dependencies
#include "mpi.h"
//#include "mkl.h"

#include <iostream>
#include <fstream>

// Internal dependencies
#include "peigen.h"

#define USE_PROFILER
#define USE_BOOST_CHRONO
#include "tic-toc-profiler.hpp"

using namespace std;
using namespace Eigen;
using namespace peigen;

// Dymode helper functions
#include "ColumnNorm.h"
//#include "H5inDMD.h"
#include "H5import.h"
#include "GEOinDMD.h"
#include "ModeSort.h"
#include "options.h"

#include <boost/math/constants/constants.hpp>

using namespace phdfp;

#endif

#include "col2ensight.h"
#include "row2text.h"



int main(int argc, char* argv[])
{
//cout << "hello! " << endl;
	//int mkl_res = mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
	MPI::Init();
//cout << "Howdy? "<<endl;
	// Deal with input parameters
	options opt(argc, argv);
//cout << "Rock n roll " << endl;
/*
	cout << "filename: " << opt.filename << endl;
	cout << "datasets: " << opt.variables[0] << endl;
	cout << "n files: " << opt.nfiles << endl;
	cout << "stride: " << opt.stride << endl;
	cout << "block: " << opt.sblock << endl;
	cout << "Eigen: " << opt.eigSolver << endl;
	cout << "res: " << opt.dispResiduals << endl;
	cout << "singulars: " << opt.nsingulars << endl;
	cout << "outdir: " << opt.outdir << endl;
	//cout << "sort: " << opt.sortMeth;
	cout << "pod: " << opt.npod << endl;
	cout << "dmd: " << opt.nmodes << endl;
*/	
//cout << "couillon "<<endl;

	int rank, numtasks;
	rank = MPI::COMM_WORLD.Get_rank();
	numtasks = MPI::COMM_WORLD.Get_size();
	bool ROOT = (rank == 0);

	DoProfiler prof;

	// Create the BLACS grid
	BLACS::init(numtasks);
	//if (BLACS::ROOT)
	//std::cout << "BLACS Initialized. MKL says " << mkl_res << endl << flush;
	//BLACS::printGrid();

	MPI::COMM_WORLD.Barrier(); // For printing purposes

	prof.tic("Dymode");

	//geofilereader georead(opt.geofile);

	/////**************************************************************************************************/
	/////*----------------------------------       READ THE DATA      ------------------------------------*/
	/////**************************************************************************************************/
	if (ROOT)
		cout << endl << "         READING DATA" << endl << "******************************" << endl;
	MPI::COMM_WORLD.Barrier(); // For printing purposes

	// Use all the processes for faster IO, regardless of them being used in the process grid
	prof.tic("Read");

	datasetreader dreader(opt.nfiles, opt.filename);

	dreader.read(opt.variables);
//cout << "SUCCESS " <<endl;
	SharedMatrix<MatrixXd> snaps(dreader.createShared(opt.sblock, opt.sblock, opt.stride));

	if (ROOT)
		prof.toc("Read", "\nReading completed in (s): ");
	else
		prof.toc("Read");

	if (ROOT)
	{
		cout << endl << "The snapshot matrix is    " << snaps.rows() << " by " << snaps.cols() << endl;
		cout << "The block size is         " << snaps.cblock() << endl;
	}

	/////**************************************************************************************************/
	/////*----------------------------------       /READ THE DATA      -----------------------------------*/
	/////**************************************************************************************************/
	// Discard inactive processes immediately, this is to avoid crashes caused e.g. by inactive process calling barrier(ctxt, All), or descinit()
	// Note: If needed, use the MPI::Intracomm BLACS::COMM_ACTIVE to avoid deadlocks with incative processes
	if (BLACS::active) // Only processes that have a place in the grid
	{
		///**************************************************************************************************/
		///*-------------------------------------     START OF DMD     -------------------------------------*/
		///**************************************************************************************************/

		prof.tic("Computations");
		/////**************************************************************************************************/
		/////*-------------------------------------       DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << endl << " SINGULAR VALUE DECOMPOSITION " << endl << "******************************" << endl;
		BLACS::COMM_ACTIVE.Barrier(); // For printing purposes

		if (ROOT)
			cout << "Calling ScaLAPACK" << endl << "=================" << endl;
		prof.tic("SVD");
		ScaSVD<MatrixXd> svd(snaps.block(0, 0, snaps.rows(), snaps.cols() - 1), true, true);
		snaps.clear();


		BLACS::COMM_ACTIVE.Barrier();
		if (ROOT)
		{
			if (opt.nsingulars > 0)
			{
				cout.precision(std::numeric_limits< double >::digits10);
				int nsings = min((int)svd.singularValues.rows(), opt.nsingulars);
				cout << "First " << nsings << " singular values: " << endl << svd.singularValues.col(0).head(nsings).transpose() << endl;
				
				cout << svd.singularValues.col(0).head(nsings).transpose() / svd.singularValues.col(0).sum() << endl << endl << flush;
				std::cout.copyfmt(std::ios(NULL));
			}
		}

		/*if (ROOT)
			cout << "HERE COMES U" << endl << flush;
			cout << svd.matrixU;*/

		/*if (ROOT)
			cout << "HERE COMES Vt" << endl << flush;
			cout << svd.matrixVt;*/

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualSVD");
			double r_svd = svd.global_residual(snaps.block(0, 0, snaps.rows(), snaps.cols() - 1));
			prof.toc("residualSVD");
			if (ROOT)
				cout << "Residual from SVD:    " << r_svd << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
		}

		if (ROOT)
			prof.toc("SVD", "\nSVD completed in (s): ");
		else
			prof.toc("SVD");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*-------------------------------------      /DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*-------------------------------      Print singular values      --------------------------------*/
		/////**************************************************************************************************/

		if (ROOT)
		{
			cout << endl << "Saving singular values...";
			prof.tic("WriteSingulars");

			std::ofstream s(opt.outdir + "singulars.txt");
			s.precision(std::numeric_limits< double >::digits10);
			if (s.is_open())
			{
				s << svd.singularValues << '\n';
				s.close();
				cout << "\tDONE." << endl;
			}
			else
			{
				cout << "\tError, could not open " << opt.outdir + "singulars.txt" << endl;
			}
			prof.toc("WriteSingulars");
		}

		/////**************************************************************************************************/
		/////*------------------------------      /Print singular values      --------------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*---------------------------       PRINT SOME POD MODES TO HDF5      ----------------------------*/
		/////**************************************************************************************************/
		if (opt.npod > 0)
		{
			prof.tic("SavePODModes");

			string podFile = opt.outdir + "POD.h5";
			stringstream variables_gold;

			//Write the time coefficients
			MatrixXd TimeSeries;
					
			if (BLACS::myrank == 0)
			{
				TimeSeries.resize(svd.matrixVt.cols(), opt.npod);
			}

			//cout << "oh oui" << endl;

			for (int m = 0; m < opt.npod; ++m)
			{
			  //cout << "m = " << m << endl;
				MatrixXd times = row2single(svd.matrixVt, m);
				//cout << "oh non" << endl;
				if (BLACS::myrank == 0)
				{
					TimeSeries.col(m) = times;

				}

			}

			if (BLACS::myrank == 0)
			{
				hid_t file = H5Fcreate(podFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
				
				hsize_t dims[2];
				dims[0] = opt.npod;
				dims[1] = svd.matrixVt.cols();
				hid_t memory_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
				hid_t file_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
//cout << TimeSeries.block(0, 0, 10, 5) << endl;
				hid_t dset = H5Dcreate( file, "/TimeSeries", H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//cout << TimeSeries << endl;
				H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, TimeSeries.data() );

				H5Dclose(dset);
				H5Fclose(file);
			}

			//Write the POD modes
			MatrixXd POD;
	//cout << "TESTING MIC" << endl;	
//BLACS::COMM_ACTIVE.Barrier(); // For debug printing purposes			
			if (BLACS::myrank == 0)
			{
				POD.resize(svd.matrixU.rows(), opt.npod);
			}
			for (int m = 0; m < opt.npod; ++m)
			{
//BLACS::COMM_ACTIVE.Barrier(); // For debug printing purposes		
//cout << "m' = " << m << endl;
				MatrixXd mode = col2single(svd.matrixU, m);
//cout << "oh noes" << endl;
				if (BLACS::myrank == 0)
				{
					POD.col(m) = mode;
				}

			}
BLACS::COMM_ACTIVE.Barrier(); // For debug printing purposes
//cout << "TESTING John" << endl;	

			if (BLACS::myrank == 0)
			{
				hid_t file = H5Fopen(podFile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);


				double Etot = svd.singularValues.sum();
				MatrixXd Energies = svd.singularValues / Etot;

				hsize_t dims[2];
				dims[0] = Energies.rows();
				dims[1] = Energies.cols();
				hid_t memory_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
				hid_t file_space = H5Screate_simple(/*rank*/ 2, dims, NULL);

				hid_t dset = H5Dcreate( file, "/Energies", H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
				H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Energies.data() );


				for (int field = 0; field < opt.variables.size(); ++field)
				{
					vector<string> dsetpaths;
					boost::split(dsetpaths, opt.variables[field], boost::is_any_of("/\\"));
					string groupname = "";
					for (int g = 1; g < (dsetpaths.size() -1); ++g)
					{
	//cout << dsetpaths[g]<< endl;
		groupname = groupname + "/" + dsetpaths[g];
						//group = H5Gcreate(group, dsetpaths[g].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					}
	//cout << groupname << endl;
					hid_t group = H5Gcreate(file, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//cout << groupname << "was created " << group << endl;
					MatrixXd buffer = POD.block(dreader.ranges[opt.variables[field]].beg, 0, 
									dreader.ranges[opt.variables[field]].end - dreader.ranges[opt.variables[field]].beg +1, POD.cols());

					
					dims[0] = buffer.cols();
					dims[1] = buffer.rows();
					memory_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
					file_space = H5Screate_simple(/*rank*/ 2, dims, NULL);

					hid_t dset = H5Dcreate( file, opt.variables[field].c_str(), H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	//cout << TimeSeries << endl;
					H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );

					H5Dclose(dset);
				}
				H5Fclose(file);
			}

			prof.toc("SavePODModes");
		}
		/////**************************************************************************************************/
		/////*----------------------------       /PRINT SOME POD MODES TO HDF5      --------------------------*/
		/////**************************************************************************************************/


		///////**************************************************************************************************/
		///////*-----------------------------      DO B := Ut * S2 * V * SIG+     ------------------------------*/
		///////**************************************************************************************************/
		if (ROOT)
			cout << endl << "      Ut * S2 * V * Sig+      " << endl << "******************************" << endl;
		BLACS::COMM_ACTIVE.Barrier(); // For printing purposes

		prof.tic("MultiplyB");

		SharedMatrix<MatrixXd> tmpMat = svd.matrixU.transpose() * snaps.block(0, 1, snaps.rows(), snaps.cols() - 1);

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
		double pinv_tol = std::numeric_limits<double>::epsilon() * snaps.rows() * svd.singularValues(0);
		if (ROOT)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			cout << "pseudo inverse tolerance: " << pinv_tol << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
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

		/*if (ROOT)
			cout << "HERE COMES Ut M V SIG+" << endl << flush;
			cout << B;*/

		if (ROOT)
			prof.toc("MultiplyB", "\nB matrix computed in (s): ");
		else
			prof.toc("MultiplyB");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*------------------------------     /DO B := Ut * M * V * SIG+     ------------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*------------------------------     DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << endl << "         EIGEN PROBLEM        " << endl << "******************************" << endl;
		BLACS::COMM_ACTIVE.Barrier(); // For printing purposes

		prof.tic("EigenProblem");

		ScaEigenSolver<MatrixXd> eig(B, true, opt.eigSolver, opt.dispResiduals);

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualEig");
			double r_eig = eig.global_residual(B);
			prof.toc("residualEig");
			if (ROOT)
				cout << "Residual from Eigen problem:        " << r_eig << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
		}

		// Matrix of eigen vectors
		SharedMatrix<MatrixXcd> X = eig.eigenVectors();
		MatrixXcd lambdas = eig.eigenValues();

		//cout << X << endl;
		if (ROOT)
			prof.toc("EigenProblem", "\nEigen problem solved in (s):        ");
		else
			prof.toc("EigenProblem");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*------------------------------    /DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/

		/////**************************************************************************************************/
		/////*----------------------------     COMPUTE THE RESCALED SPECTRUM     -----------------------------*/
		/////**************************************************************************************************/

		MatrixXd ScaledAmplitudes;
		if (opt.sortMeth.stype == scaled)
		{
			if (ROOT)
			{	cout << "Computing the rescaled spectrum...             " << endl;
			}

			SharedMatrix<MatrixXcd> ScaledModes = svd.matrixVt.cast<std::complex<double> >().transpose();
			ScaledModes.local_matrix = ScaledModes.local_matrix * SIGplus;

			SharedMatrix<MatrixXcd> tmp = ScaledModes * X;

			MatrixXd ScalingValues = ColumnNorm(tmp);
			tmp.ColScale(ScalingValues.cast<std::complex<double> >().cwiseInverse());
			ScaledModes = (snaps.cast<std::complex<double> >().block(0, 0, snaps.rows(), snaps.cols() - 1)) * tmp;

			ScaledAmplitudes = ColumnNorm(ScaledModes);
		}
		
		/////**************************************************************************************************/
		/////*----------------------------     /COMPUTE THE RESCALED SPECTRUM     ----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*------------------------------      DO A LINEAR SYSTEM SOLVE      ------------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << endl << "         LINEAR SYSTEM        " << endl << "******************************" << endl;
		BLACS::COMM_ACTIVE.Barrier(); // For printing purposes

		prof.tic("LinearSolve");
		//cout << "(" << BLACS::myrank << ")" << endl;
		prof.tic("FormRHS");

		if (ROOT)
			cout << "Preparing the right-hand-side...";

		SharedMatrix<MatrixXd> rhs = svd.matrixU.transpose() * snaps.block(0, 0/*Nt - 1*/, snaps.rows(), 1);
		svd.matrixU.clear();
		snaps.clear();

		SharedMatrix<MatrixXcd> rhsZ = rhs.cast<std::complex<double> >();
		prof.toc("FormRHS");

		if (ROOT)
			cout << "\tDONE" << endl;

		if (ROOT)
			cout << "Preparing the system...         ";

		BLACS::COMM_ACTIVE.Barrier();

		prof.tic("FormSystem");
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
		prof.toc("FormSystem");

		BLACS::COMM_ACTIVE.Barrier();

		if (ROOT)
			cout << "\tDONE" << endl << endl;

		BLACS::COMM_ACTIVE.Barrier();

		//ScaSolve<MatrixXcd> solver(X, rhsZ, peigen::EigenSVD);

		if (ROOT)
			cout << "Calling ScaLAPACK" << endl << "=================" << endl;


		BLACS::COMM_ACTIVE.Barrier();


		prof.tic("SolveSystem");
		ScaSolve<MatrixXd> solver(System, rhs, peigen::pxgesvx);
		prof.toc("SolveSystem");

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualLin");
			double r_lin = solver.residual(System, rhs);
			prof.toc("residualLin");
			if (ROOT)
				cout << "Residual from the linear system: \t" << r_lin << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
		}

		//cout << BLACS::myrank << ", lambdas: " << lambdas << endl << flush;
		//cout << BLACS::myrank << ", solution: " << solver.solution.local_matrix << endl << flush;

		if (ROOT)
			cout << "Reconstructing the weights...     ";

		prof.tic("FormWeights");
		// Reconstitute the solution to the original system
		SharedMatrix<MatrixXcd> weights(solver.solution.rows(), solver.solution.cols(), solver.solution.rblock(), solver.solution.cblock());
		if (weights.local_matrix.cols() > 0)
		{
			for (int k = 0; k < lambdas.rows(); ++k)
			{
				// First item of a conjugate pair
				if (lambdas(k, 0).imag() != 0)
				{
					if (BLACS::indxg2p(k, weights.rblock(), BLACS::grid_rows) == BLACS::myrow)
					{
						int l = BLACS::indxg2l(k, weights.rblock(), BLACS::grid_rows);
						if (BLACS::indxg2p(k + 1, weights.rblock(), BLACS::grid_rows) == BLACS::myrow)
						{
							//cout << BLACS::myrank << ", I have a pair " << k << endl << flush;
							weights.local_matrix.row(l).real() = solver.solution.local_matrix.row(l);
							weights.local_matrix.row(l).imag() = solver.solution.local_matrix.row(l + 1);
							weights.local_matrix.row(l + 1) = weights.local_matrix.row(l).conjugate();
						}
						else // have the 1st one but not the 2nd one
						{
							int ownernext = BLACS::indxg2p(k + 1, weights.rblock(), BLACS::grid_rows);
							ownernext = BLACS::Cblacs_pnum(BLACS::ctxt, ownernext, BLACS::mycol);
							MatrixXd re = solver.solution.local_matrix.row(l);
							MatrixXd im(1, solver.solution.local_matrix.cols());
							//cout << BLACS::myrank << ", sending " << k << " to " << ownernext << " with tag " << ownernext << endl << flush;
							BLACS::COMM_ACTIVE.Send(re.data(), re.cols(), MPI::DOUBLE, ownernext, ownernext);
							BLACS::COMM_ACTIVE.Recv(im.data(), im.cols(), MPI::DOUBLE, ownernext, ownernext);

							weights.local_matrix.row(l).real() = solver.solution.local_matrix.row(l);
							weights.local_matrix.row(l).imag() = im;
						}
					}
					else
					{
						if (BLACS::indxg2p(k + 1, weights.rblock(), BLACS::grid_rows) == BLACS::myrow)
						{
							int l = BLACS::indxg2l(k + 1, weights.rblock(), BLACS::grid_rows);

							int ownerprev = BLACS::indxg2p(k, weights.rblock(), BLACS::grid_rows);
							ownerprev = BLACS::Cblacs_pnum(BLACS::ctxt, ownerprev, BLACS::mycol);
							MatrixXd re(1, solver.solution.local_matrix.cols());
							MatrixXd im = solver.solution.local_matrix.row(l);
							//cout << BLACS::myrank << ", receiving " << k << " from " << ownerprev << " with tag " << BLACS::myrank << endl << flush;
							BLACS::COMM_ACTIVE.Recv(re.data(), re.cols(), MPI::DOUBLE, ownerprev, BLACS::myrank);
							BLACS::COMM_ACTIVE.Send(im.data(), im.cols(), MPI::DOUBLE, ownerprev, BLACS::myrank);

							weights.local_matrix.row(l).real() = re;
							weights.local_matrix.row(l).imag() = -solver.solution.local_matrix.row(l);
						}
					}
					++k;
				}
				else // real eigenvalue
				{
					if (BLACS::indxg2p(k, weights.rblock(), BLACS::grid_rows) == BLACS::myrow)
					{
						//cout << BLACS::myrank << ", I have a real one " << k << endl << flush;
						int l = BLACS::indxg2l(k, weights.rblock(), BLACS::grid_rows);
						weights.local_matrix.row(l) = solver.solution.local_matrix.row(l).cast<complex<double>>();
					}
				}
			}
		}
		prof.toc("FormWeights");

		if (ROOT)
			cout << "\tDONE" << endl;


		/////**************************************************************************************************/
		/////*------------------------------      /DO A LINEAR SYSTEM SOLVE      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*-----------------------------      APPLY WEIGHT TO THE MODES      ------------------------------*/
		/////**************************************************************************************************/

		if (ROOT)
			cout << "Creating the modes...             ";

		SharedMatrix<MatrixXcd> Modes = svd.matrixU.cast<std::complex<double> >() * X;

		if (ROOT)
			cout << "\tDONE" << endl;

		if (ROOT)
			cout << "Scaling the modes...              ";
		Modes.ColScale(weights);
		if (ROOT)
			cout << "\tDONE" << endl;

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualLin");

			SharedMatrix<MatrixXcd> Vandermonde = vander<MatrixXcd>(lambdas, Modes.cols() + 1, Modes.rblock(), Modes.cblock());
			//cout << Vandermonde << endl;

			SharedMatrix<MatrixXcd> reconstruct = snaps.cast<complex<double>>();
			reconstruct.pgemm(1., Modes, Vandermonde, -1.);
			//cout << reconstruct << endl;

			reconstruct.block(0, 0, snaps.rows(), snaps.cols() - 1);
			double r_loc = (reconstruct.localBlock().rows() * reconstruct.localBlock().cols()) > 0
				? reconstruct.localBlock().cwiseAbs().maxCoeff()
				: -1;
			double r;
			BLACS::COMM_ACTIVE.Reduce(&r_loc, &r, 1, MPI::DOUBLE, MPI::MAX, 0);

			if (ROOT)
				cout << "Residual from Modes:                \t" << r << endl;

			reconstruct.block(0, snaps.cols() - 1, snaps.rows(), 1);
			r_loc = (reconstruct.localBlock().rows() * reconstruct.localBlock().cols()) > 0
				? reconstruct.localBlock().cwiseAbs().maxCoeff()
				: -1;

			BLACS::COMM_ACTIVE.Reduce(&r_loc, &r, 1, MPI::DOUBLE, MPI::MAX, 0);
			if (ROOT)
				cout << "Residual from Modes (last snapshot): \t" << r << endl;
			prof.toc("residualLin");


			std::cout.copyfmt(std::ios(NULL));
		}

		if (ROOT)
			prof.toc("LinearSolve", "\nLinear system solved in (s):    ");
		else
			prof.toc("LinearSolve");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*-----------------------------      /APPLY WEIGHT TO THE MODES      -----------------------------*/
		/////**************************************************************************************************/


		if (ROOT)
			cout << endl << "          SAVING DATA" << endl << "******************************" << endl;
		/////**************************************************************************************************/
		/////*------------------------------      Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << "Computing the modes' norm...";
		prof.tic("Energy");
		MatrixXd amplitudes = ColumnNorm(Modes);
		prof.toc("Energy");
		if (ROOT)
			cout << "\tDONE";

		/////**************************************************************************************************/
		/////*-----------------------------      /Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/
		prof.toc("Computations");

		/////**************************************************************************************************/
		/////*----------------------------------      Print light data      ----------------------------------*/
		/////**************************************************************************************************/

		if (ROOT)
		{
			prof.tic("WriteLight");
			cout << endl << "Saving spectrum...          ";

			Matrix<double, Dynamic, 2, RowMajor> spectrum(snaps.cols() - 1, 2);
			spectrum.col(0) = lambdas.imag().binaryExpr(lambdas.real(), std::ptr_fun(atan2<double, double>));
			spectrum.col(1) = amplitudes.transpose();

			std::ofstream s(opt.outdir + "spectrum.txt");
			s.precision(std::numeric_limits< double >::digits10);
			if (s.is_open())
			{
				s << spectrum << '\n';
				s.close();
				cout << "\tDONE" << endl;
			}
			else
			{
				cout << "\tError, could not open " << opt.outdir + "spectrum.txt" << endl;
			}

			cout << "Saving eigenvalues...       ";
			std::ofstream l(opt.outdir + "eigenvalues.txt");
			l.precision(std::numeric_limits< double >::digits10);
			if (l.is_open())
			{
				l << lambdas << '\n';
				l.close();
				cout << "\tDONE" << endl;
			}
			else
			{
				cout << "\tError, could not open " << opt.outdir + "eigenvalues.txt" << endl;
			}


			if (opt.sortMeth.stype == scaled)
			{
				cout << "Saving scaled spectrum...       ";
				std::ofstream scaledfile(opt.outdir + "scaled_spectrum.txt");
				scaledfile.precision(std::numeric_limits< double >::digits10);
				if (scaledfile.is_open())
				{
					scaledfile << ScaledAmplitudes << '\n';
					scaledfile.close();
					cout << "\tDONE" << endl;
				}
				else
				{
					cout << "\tError, could not open " << opt.outdir + "scaled_spectrum.txt" << endl;
				}
			}

			prof.toc("WriteLight", "\nLight data saved in (s): ");
			cout << endl;
		}

		/////**************************************************************************************************/
		/////*----------------------------------      /Print light data      ---------------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*---------------------------       PRINT SOME DMD MODES TO HDF5      ----------------------------*/
		/////**************************************************************************************************/
		prof.tic("SaveModes");
		BLACS::COMM_ACTIVE.Barrier();
		// Find which column of Modes has the most energy
		MatrixXd::Index i_mode, x_mode;
		stringstream variables_gold;

		MatrixXd SortingAmplitude;
		if (opt.sortMeth.stype == scaled)
		{
			SortingAmplitude = ScaledAmplitudes;
		}
		else
		{
			SortingAmplitude = amplitudes;
		}

		//cout << BLACS::myrank << " amplitudes " << amplitudes << endl << endl;
		ModeSort<MatrixXcd> sorted(Modes, lambdas, SortingAmplitude, svd.singularValues, opt.sortMeth, opt.nmodes);
		const MatrixXi indices = sorted.orderedIdx;

		MatrixXcd DMD;
		if (ROOT)
		{
			DMD.resize(Modes.rows(), indices.cols());
cout << "DMD = " << DMD.rows()<< "x" << DMD.cols() << endl;
		}


		MatrixXd mean_energy(1, amplitudes.cols());
		double Etot = 0;
		for (int m = 0; m < amplitudes.cols(); ++m)
		{
			if (abs(lambdas(m, 0)) == 1)
				mean_energy(0, m) = amplitudes(0, m);
			else
				mean_energy(0, m) = amplitudes(0, m)*(1 - std::pow(abs(lambdas(m, 0)), lambdas.rows())) / (lambdas.rows() * (1 - abs(lambdas(m, 0))));
			Etot += mean_energy(0, m);
		}
		

		MatrixXd Energies(1, indices.cols());
		MatrixXd Frequencies(1, indices.cols());
		for (int m = 0; m < indices.cols(); ++m)
		{
//cout << "m = " << m << endl;
			const int i_mode = indices(0, m);	
//cout << "i_mode = " << i_mode << endl;
			Energies(0, m) = mean_energy(0, i_mode) / Etot;
			Frequencies(0, m) = atan2(lambdas(i_mode).imag(), lambdas(i_mode).real()) / (2*opt.tstep*boost::math::constants::pi<double>());

			MatrixXcd curr_mode = col2single(Modes, i_mode);
			if (ROOT)
			{
//cout << "yeapeeeh" << endl;
//cout << "oh yeahs: " << Modes.local_matrix.block(dreader.ranges[opt.variables[1]].beg,i_mode, 10, 1) << endl;
//cout << "oh noes: "  << endl;
				DMD.col(m) = curr_mode;
			}		
		}
//cout << "lambdas: "  << endl<< lambdas.cwiseAbs() <<endl;
		if (ROOT)
		{
			string dmdFile = opt.outdir + "DMD.h5";
			hid_t file = H5Fcreate(dmdFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


			hsize_t dims[2];
			dims[0] = Energies.rows();
			dims[1] = Energies.cols();
			hid_t memory_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
			hid_t file_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
			hid_t dset = H5Dcreate( file, "/Energies", H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
			H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Energies.data() );

			dset = H5Dcreate( file, "/Frequencies", H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
			H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Frequencies.data() );

			for (int field = 0; field < opt.variables.size(); ++field)
			{
				vector<string> dsetpaths;
				boost::split(dsetpaths, opt.variables[field], boost::is_any_of("/\\"));
				string groupname = "";
				for (int g = 1; g < (dsetpaths.size() -1); ++g)
				{
					groupname = groupname + "/" + dsetpaths[g];
				}
				hid_t group = H5Gcreate(file, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//cout << "taking a block of " << dreader.ranges[opt.variables[field]].beg << ",0; " << dreader.ranges[opt.variables[field]].end - dreader.ranges[opt.variables[field]].beg +1 << "," <<DMD.cols() << endl;

	//cout << DMD.block(dreader.ranges[opt.variables[field]].beg, 0, 10, 2) << endl << endl;
	//cout << Modes.local_matrix.block(dreader.ranges[opt.variables[field]].beg, indices(0,0), 10, 1) << endl << endl;
	//cout << Modes.local_matrix.block(dreader.ranges[opt.variables[field]].beg, indices(0,1), 10, 1) << endl << endl;

				MatrixXd buffer = DMD.block(dreader.ranges[opt.variables[field]].beg, 0, 
								dreader.ranges[opt.variables[field]].end - dreader.ranges[opt.variables[field]].beg +1, DMD.cols()).cwiseAbs().cast<double>();

				dims[0] = buffer.cols();
				dims[1] = buffer.rows();
				memory_space = H5Screate_simple(/*rank*/ 2, dims, NULL);
				file_space = H5Screate_simple(/*rank*/ 2, dims, NULL);

				string dsetname = opt.variables[field] + "_Magnitude";

//cout << "mag " << endl << buffer.block(0, 0, 10, 2) << endl << endl;

				dset = H5Dcreate( file, dsetname.c_str(), H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	//cout << TimeSeries << endl;
				H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );

				buffer = DMD.block(dreader.ranges[opt.variables[field]].beg, 0, 
								dreader.ranges[opt.variables[field]].end - dreader.ranges[opt.variables[field]].beg +1, DMD.cols()).imag().binaryExpr(DMD.block(dreader.ranges[opt.variables[field]].beg, 0, 
								dreader.ranges[opt.variables[field]].end - dreader.ranges[opt.variables[field]].beg +1, DMD.cols()).real(), std::ptr_fun(atan2<double, double>)).cast<double>();

//cout << "phase " << endl << buffer.block(0, 0, 10, 2) << endl << endl;
				
				H5Dclose(dset);

				dsetname = opt.variables[field] + "_Phase";
				dset = H5Dcreate( file, dsetname.c_str(), H5T_IEEE_F64LE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
				H5Dwrite(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data() );
				H5Dclose(dset);

				
			}
			H5Fclose(file);
		}

		////////////////////////
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		if (ROOT)
			prof.toc("SaveModes", "\nModes saved in (s):      ");
		else
			prof.toc("SaveModes");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*------------------------------      /PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			prof.toc("Dymode", "Dymode completed in (s): ");
		stringstream profile_data;
		profile_data << opt.outdir << "profiler" << "-" << rank << ".yml";

		//if (BLACS::myrank == 0)
		prof.dump(profile_data.str());

		BLACS::COMM_ACTIVE.Barrier();
	} // end IF ACTIVE
	else
	{
	}

	MPI::COMM_WORLD.Barrier();

	if (ROOT)	// not having if(root) prevents segfault at the begining of the program (!?!)
		cout << endl << endl << "DYMODE OUT!" << endl;

	BLACS::finalize();
	MPI::Finalize();

	return 0;
}

