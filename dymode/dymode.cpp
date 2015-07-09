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
#include "H5inDMD.h"
#include "GEOinDMD.h"
#include "ModeSort.h"
#include "options.h"

using namespace phdfp;

#endif



void set80line(string& s)
{
	s.resize(80, ' ');
	s.back() = '\n';
}
void gold_print_header(int mode, string complex_part, string var, FILE *pFile)
{
	stringstream stext;
	string text;

	stext << complex_part << " of Mode "
		<< setfill('0') << setw(6) << mode
		<< " for " << var;
	text = stext.str();
	set80line(text);

	fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
	stext.clear();//clear any bits set
	stext.str(std::string());
}

void gold_print_values(MatrixXf values, geofilereader geo, FILE *pFile)
{
	stringstream stext;
	string text;

	int offset = 0;

	for (auto it_part = geo.parts.begin(); it_part != geo.parts.end(); ++it_part)
	{
		stext << "part";
		text = stext.str();
		set80line(text);

		fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
		stext.clear();//clear any bits set
		stext.str(std::string());

		int part_number = (*it_part).number;
		fwrite(&part_number, 1, 1 * sizeof(int), pFile);

		for (unsigned int k = 0; k < (*it_part).telements.size(); ++k)
		{
			stext << (*it_part).telements[k];
			text = stext.str();
			set80line(text);

			fwrite(text.c_str(), 1, 80 * sizeof(char), pFile);
			stext.clear();//clear any bits set
			stext.str(std::string());

			int nelems = (*it_part).nelements[k];
			if (nelems == -1)
				nelems = values.rows();
			fwrite(values.data() + offset, sizeof(float), nelems, pFile);
			offset += (*it_part).nelements[k];
		}
	}
}


int main(int argc, char* argv[])
{
	//int mkl_res = mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
	MPI::Init();

	// Deal with input parameters
	options opt(argc, argv);

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

	geofilereader georead(opt.geofile);

	/////**************************************************************************************************/
	/////*----------------------------------       READ THE DATA      ------------------------------------*/
	/////**************************************************************************************************/
	if (ROOT)
		cout << endl << "         READING DATA" << endl << "******************************" << endl;
	MPI::COMM_WORLD.Barrier(); // For printing purposes

	// Use all the processes for faster IO, regardless of them being used in the process grid
	prof.tic("Read");

	datasetreader dreader(opt.nfiles, opt.filename);

	dreader.read(opt.dataset, opt.variables);

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
				cout << "First " << nsings << " singular values: " << endl << svd.singularValues.col(0).head(nsings).transpose() << endl << endl << flush;
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

			SharedMatrix<MatrixXcd> ScaledModes = svd.matrixVt.transpose();
			ScaledModes = ScaledModes.local_matrix * SIGplus;
			ScaledModes = ScaledModes * X;

			MatrixXd ScalingValues = ColumnNorm(ScaledModes);
			ScaledModes.ColScale(ScalingValues.cwiseInverse);
			ScaledModes = snaps.block(0, 0, snaps.rows(), snaps.cols() - 1) * ScaledModes;

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
		/////*------------------------------       PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/
		prof.tic("SaveModes");
		BLACS::COMM_ACTIVE.Barrier();
		// Find which column of Modes has the most energy
		MatrixXd::Index i_mode, x_mode;
		stringstream variables_gold;

		MatrixXd SortingAmplitude:
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
				if (BLACS::mycol == BLACS::indxg2p(i_mode, Modes.cblock(), BLACS::grid_cols))
				{
					int i_loc = BLACS::indxg2l(i_mode, Modes.cblock(), BLACS::grid_cols);
					MatrixXf exportdata;


					// Gather the global column on row 0
					Matrix<MPI::Request, Dynamic, Dynamic> Irecv_requests;
					Matrix<Matrix<float, Dynamic, Dynamic>, Dynamic, 1> RecvBuffer;
					if (BLACS::myrow == 0)
					{
						exportdata.resize(Modes.rows(), 2);
						Irecv_requests.resize(BLACS::grid_rows, 1);
						RecvBuffer.resize(BLACS::grid_rows, 1);

						for (int r = 0; r < BLACS::grid_rows; ++r)
						{
							int size = BLACS::peigen_numroc(Modes.rows(), Modes.rblock(), r, 0, BLACS::grid_rows);
							RecvBuffer(r, 0).resize(size, 2);
							int powner = BLACS::Cblacs_pnum(BLACS::ctxt, r, BLACS::mycol);
							Irecv_requests(r, 0) = BLACS::COMM_ACTIVE.Irecv(RecvBuffer(r, 0).data(), size * 2, MPI::FLOAT, powner, 1/*tag*/);
						}
					}

					Matrix<float, Dynamic, 2> SendBuffer(Modes.local_matrix.rows(), 2);
					SendBuffer.col(0) = Modes.local_matrix.col(i_loc).cwiseAbs().cast<float>();
					SendBuffer.col(1) = Modes.local_matrix.col(i_loc).imag().binaryExpr(Modes.local_matrix.col(i_loc).real(), std::ptr_fun(atan2<double, double>)).cast<float>();
					int col_root = BLACS::Cblacs_pnum(BLACS::ctxt, 0, BLACS::mycol);
					BLACS::COMM_ACTIVE.Send(SendBuffer.data(), SendBuffer.rows() * SendBuffer.cols(), MPI::FLOAT, col_root, 1 /*tag*/);

					if (BLACS::myrow == 0)
					{
						MPI::Request::Waitall(BLACS::grid_rows, Irecv_requests.data());

						// Combine the buffers
						for (int rb = 0; rb < ceil((double)exportdata.rows() / Modes.rblock()); rb++)
						{
							int roffset = Modes.rblock() * floor(rb / BLACS::grid_rows);
							int _nrows = min(Modes.rblock(), (int)(exportdata.rows() - rb*Modes.rblock()));
							int pr_owner = rb % BLACS::grid_rows;
							exportdata.block(rb*Modes.rblock(), 0, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 0, _nrows, 1);
							exportdata.block(rb*Modes.rblock(), 1, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 1, _nrows, 1);
						}
					}

					// Print to disk from row 0
					if (BLACS::myrow == 0)
					{
						cout << "(" << BLACS::myrank << ") " << " writing mode " << m << "...";
						FILE *pFile;

						int offset_gold = 0;
						MatrixXf values;
						for (string var : opt.variables)
						{
							if (!(var == "null"))
							{
								stringstream filenameRE;
								filenameRE << "mode" << setfill('0') << setw(6) << m << "." << var << ".abs";

								pFile = fopen((opt.outdir + filenameRE.str()).c_str(), "wb");

								gold_print_header(m, "Module", var, pFile);

								values = exportdata.block(offset_gold, 0, dreader.Np, 1);
								gold_print_values(values, georead, pFile);

								fclose(pFile);

								stringstream filenameIM;
								filenameIM << "mode" << setfill('0') << setw(6) << m << "." << var << ".ang";

								pFile = fopen((opt.outdir + filenameIM.str()).c_str(), "wb");

								gold_print_header(m, "Angle", var, pFile);

								values = exportdata.block(offset_gold, 1, dreader.Np, 1);
								gold_print_values(values, georead, pFile);

								fclose(pFile);

								offset_gold += dreader.Np;
							}
						}
						cout << "\tDONE" << endl;
					}
				}

				// Add all modes/variables to the list of variables
				if (BLACS::myrank == 0)
				{
					for (string var : opt.variables)
					{
						if (!(var == "null"))
						{
							variables_gold << "scalar per element: " << var << m << "abs " << "mode" << setfill('0') << setw(6) << m << "." << var << ".abs" << endl;
							variables_gold << "scalar per element: " << var << m << "ang " << "mode" << setfill('0') << setw(6) << m << "." << var << ".ang" << endl;
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

