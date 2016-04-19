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
			stringstream variables_gold;
			MatrixXd TimeSeries;

			if (BLACS::myrank == 0)
			{
				TimeSeries.resize(svd.matrixVt.cols(), opt.npod);
			}

			for (int m = 0; m < opt.npod; ++m)
			{

				BLACS::COMM_ACTIVE.Barrier();

				stringstream modename;
				modename << "pod" << setfill('0') << setw(6) << m;
				//col2ensight(svd.matrixU, m, modename.str(), false, georead, dreader, opt);

				// Add all modes/variables to the list of variables
				if (BLACS::myrank == 0)
				{
					for (string var : opt.variables)
					{
						if (!(var == "null"))
						{
							variables_gold << "scalar per element: " << var << m << " " << modename.str() << "." << var << endl;
						}
					}
				}

				MatrixXd times = row2single(svd.matrixVt, m);
				if (BLACS::myrank == 0)
				{
					TimeSeries.col(m) = times;
				}

			}

			// Write the .case file
			if (BLACS::myrank == 0)
			{
				vector<string> geopath;
				boost::split(geopath, opt.geofile, boost::is_any_of("/\\"));

				std::ofstream ofs(opt.outdir + "pod.case", std::ofstream::out);
				ofs << "FORMAT" << endl
					<< "type: ensight gold" << endl
					<< "GEOMETRY" << endl
					<< "model: " << geopath.back() << endl
					<< "VARIABLE" << endl
					<< variables_gold.str()
					<< "TIME" << endl
					<< "time set: 1 \nnumber of steps: 1 \nfilename start number: 0 \nfilename increment: 1 \ntime values: \n0" << endl;
				ofs.close();
			}

			// Write the time coefficients
			if (BLACS::myrank == 0)
			{
				std::ofstream ofs(opt.outdir + "time_coeffs.txt", std::ofstream::out);
				ofs.precision(std::numeric_limits< double >::digits10);
				ofs << TimeSeries;
				ofs.close();
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

		for (int m = 0; m < indices.cols(); ++m)
		{
			const int i_mode = indices(0, m);

			if (opt.sortMeth.conjugates == false && lambdas(i_mode, 0).imag() < 0)
			{	// There is no more mode in the top half plane
				m = indices.cols();
			}
			else
			{
				stringstream modename;
				modename << "dmd" << setfill('0') << setw(6) << m;
				//col2ensight(Modes, i_mode, modename.str(), true, georead, dreader, opt);
				
				// Add all modes/variables to the list of variables
				if (BLACS::myrank == 0)
				{
					for (string var : opt.variables)
					{
						if (!(var == "null"))
						{
							variables_gold << "scalar per element: " << var << m << "abs " << modename.str() << "." << var << ".abs" << endl;
							variables_gold << "scalar per element: " << var << m << "ang " << modename.str() << "." << var << ".ang" << endl;
						}
					}
				}
			}
		}

		// Write the .case file
		if (BLACS::myrank == 0)
		{
			vector<string> geopath;
			boost::split(geopath, opt.geofile, boost::is_any_of("/\\"));

			std::ofstream ofs(opt.outdir + "dmd.case", std::ofstream::out);
			ofs << "FORMAT" << endl
				<< "type: ensight gold" << endl
				<< "GEOMETRY" << endl
				<< "model: " << geopath.back() << endl
				<< "VARIABLE" << endl
				<< variables_gold.str()
				<< "TIME" << endl
				<< "time set: 1 \nnumber of steps: 1 \nfilename start number: 0 \nfilename increment: 1 \ntime values: \n0" << endl;
			ofs.close();
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

