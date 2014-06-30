// DMD.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"

#ifndef USE_PRECOMPILED_HEADER

// External dependencies
#include "mpi.h"
#include "mkl.h"

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
#include "ColumnSquaredNorm.h"
#include "H5inDMD.h"
#include "GEOinDMD.h"
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

			fwrite(values.data() + offset, sizeof(float), (*it_part).nelements[k], pFile);
			offset += (*it_part).nelements[k];
		}
	}
}


int main(int argc, char* argv[])
{
	int mkl_res = mkl_cbwr_set(MKL_CBWR_COMPATIBLE);
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
	if (BLACS::ROOT)
		std::cout << "BLACS Initialized. MKL says " << mkl_res << endl << flush;
	//BLACS::printGrid();

	MPI::COMM_WORLD.Barrier(); // For printing purposes
	


	geofilereader georead(opt.geofile);

	/////**************************************************************************************************/
	/////*----------------------------------       READ THE DATA      ------------------------------------*/
	/////**************************************************************************************************/
	if (ROOT)
		cout << endl << "****** Reading data..." << endl << endl;
	MPI::COMM_WORLD.Barrier(); // For printing purposes

	// Use all the processes for faster IO, regardless of them being used in the process grid
	prof.tic("Read");

	datasetreader dreader(opt.nfiles, opt.filename);

	dreader.read(opt.dataset, opt.variables);

	SharedMatrix<MatrixXd> snaps(dreader.createShared(6, 6, opt.stride));

	//cout << snaps << endl;

	if (ROOT)
		prof.toc("Read", "\nReading completed in (s): ");
	else
		prof.toc("Read");
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

		
		/////**************************************************************************************************/
		/////*-------------------------------------       DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << endl << "****** Computing Singular Value Decomposition..." << endl << endl;
		MPI::COMM_WORLD.Barrier(); // For printing purposes

		prof.tic("SVD");
		ScaSVD<MatrixXd> svd(snaps.block(0, 0, snaps.rows(), snaps.cols() - 1), true, true);
		snaps.clear();


		BLACS::COMM_ACTIVE.Barrier();
		if (ROOT)
		{
			if (opt.nsingulars > 0)
			{
				cout.precision(std::numeric_limits< double >::digits10);
				cout << "First " << opt.nsingulars << " singular values: " << endl << svd.singularValues.col(0).head(opt.nsingulars).transpose() << endl << endl << flush;
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
				cout << "Residual from SVD: " << r_svd << endl << flush;
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

		


		///////**************************************************************************************************/
		///////*-----------------------------      DO B := Ut * S2 * V * SIG+     ------------------------------*/
		///////**************************************************************************************************/
		if (ROOT)
			cout << endl << "****** Forming Ut * S2 * V * Sig+..." << endl << endl;
		MPI::COMM_WORLD.Barrier(); // For printing purposes
		
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
			cout << "pinv_tol: " << pinv_tol << endl << flush;
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
			cout << endl << "****** Solving the eigen problem..." << endl << endl;
		MPI::COMM_WORLD.Barrier(); // For printing purposes
		
		prof.tic("EigenProblem");

		ScaEigenSolver<MatrixXd> eig(B, true, opt.eigSolver);

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualEig");
			double r_eig = eig.global_residual(B);
			prof.toc("residualEig");
			if (ROOT)
				cout << "Residual from Eigen problem: " << r_eig << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
		}
		
		// Matrix of eigen vectors
		SharedMatrix<MatrixXcd> X = eig.eigenVectors();
		MatrixXcd lambdas = eig.eigenValues();

		if (ROOT)
			prof.toc("EigenProblem", "\nEigen problem solved in (s): ");
		else
			prof.toc("EigenProblem");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*------------------------------    /DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*------------------------------      DO A LINEAR SYSTEM SOLVE      ------------------------------*/
		/////**************************************************************************************************/
		if (ROOT)
			cout << endl << "****** Solving the linear system..." << endl << endl;
		MPI::COMM_WORLD.Barrier(); // For printing purposes
		
		prof.tic("LinearSolve");
		//cout << "(" << BLACS::myrank << ")" << endl;
		prof.tic("FormRHS");
		SharedMatrix<MatrixXd> rhs = svd.matrixU.transpose() * snaps.block(0, 0/*Nt - 1*/, snaps.rows(), 1);
		svd.matrixU.clear();
		snaps.clear();
		SharedMatrix<MatrixXcd> rhsZ = rhs.cast<std::complex<double> >();
		prof.toc("FormRHS");

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

		//ScaSolve<MatrixXcd> solver(X, rhsZ, peigen::EigenSVD);

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
				cout << "Residual from the linear system: " << r_lin << endl << flush;
			std::cout.copyfmt(std::ios(NULL));
		}

		//cout << BLACS::myrank << ", lambdas: " << lambdas << endl << flush;
		//cout << BLACS::myrank << ", solution: " << solver.solution.local_matrix << endl << flush;

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
							MatrixXd im = solver.solution.local_matrix.row(l + 1);
							//cout << BLACS::myrank << ", receiving " << k << " from " << ownerprev << " with tag " << BLACS::myrank << endl << flush;
							BLACS::COMM_ACTIVE.Recv(re.data(), re.cols(), MPI::DOUBLE, ownerprev, BLACS::myrank);
							BLACS::COMM_ACTIVE.Send(im.data(), im.cols(), MPI::DOUBLE, ownerprev, BLACS::myrank);

							weights.local_matrix.row(l + 1).real() = re;
							weights.local_matrix.row(l + 1).imag() = solver.solution.local_matrix.row(l + 1);
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

		/*if (ROOT)
		{
			cout << "HERE COME the solution" << endl << flush;
		}
		cout << weights << endl;*/



		/////**************************************************************************************************/
		/////*------------------------------      /DO A LINEAR SYSTEM SOLVE      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*-----------------------------      APPLY WEIGHT TO THE MODES      ------------------------------*/
		/////**************************************************************************************************/

		SharedMatrix<MatrixXcd> Modes = svd.matrixU.cast<std::complex<double> >() * X;

		Modes.ColScale(weights);

		if (opt.dispResiduals)
		{
			cout.precision(std::numeric_limits< double >::digits10);
			prof.tic("residualLin");

			SharedMatrix<MatrixXcd> Vandermonde = vander<MatrixXcd>(lambdas, Modes.cols(), Modes.rblock(), Modes.cblock());
			//cout << Vandermonde << endl;

			SharedMatrix<MatrixXcd> reconstruct = snaps.cast<complex<double>>().block(0, 0, snaps.rows(), snaps.cols() - 1);
			reconstruct.pgemm(1., Modes, Vandermonde, -1.);
			//cout << reconstruct << endl;

			double r_loc = reconstruct.localBlock().cwiseAbs().maxCoeff();
			double r;
			BLACS::COMM_ACTIVE.Reduce(&r_loc, &r, 1, MPI::DOUBLE, MPI::MAX, 0);
			prof.toc("residualLin");

			if (ROOT)
				cout << "Residual from Modes: " << r << endl;
			std::cout.copyfmt(std::ios(NULL));
		}
		
		if (ROOT)
			prof.toc("LinearSolve", "\nLinear system solved in (s): ");
		else
			prof.toc("LinearSolve");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*-----------------------------      /APPLY WEIGHT TO THE MODES      -----------------------------*/
		/////**************************************************************************************************/

		

		/////**************************************************************************************************/
		/////*------------------------------      Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/		
		
		prof.tic("Energy");
		MatrixXd amplitudes = ColumnSquaredNorm(Modes);
		prof.toc("Energy");


		/////**************************************************************************************************/
		/////*-----------------------------      /Compute the mode's energy      -----------------------------*/
		/////**************************************************************************************************/


		/////**************************************************************************************************/
		/////*----------------------------------      Print light data      ----------------------------------*/
		/////**************************************************************************************************/

		if (ROOT)
		{
			prof.tic("WriteLight");
			cout << endl << "****** Saving spectrum...";
			cout << endl << "*************************" << endl << " SPECTRUM IS BUGGED FOR SOME REASON. INVESTIGATE!" << endl << "*************************" << endl;
			Matrix<double, Dynamic, 2, RowMajor> spectrum(snaps.cols() - 1, 2);
			spectrum.col(0) = lambdas.imag().binaryExpr(lambdas.real(), std::ptr_fun(atan2<double, double>));
			spectrum.col(1) = amplitudes.transpose();

			std::ofstream s(opt.outdir + "spectrum.txt");
			//s.precision(std::numeric_limits< double >::digits10);
			if (s.is_open())
			{
				s << spectrum << '\n';
				s.close();
			}
			cout << "\tDONE." << endl;

			cout << endl << "****** Saving eigenvalues..." ;
			std::ofstream l(opt.outdir + "eigenvalues.txt");
			//l.precision(std::numeric_limits< double >::digits10);
			if (l.is_open())
			{
				l << lambdas << '\n';
				l.close();
			}
			cout << "\tDONE." << endl;
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

		//cout << BLACS::myrank << " amplitudes " << amplitudes << endl << endl;

		for (int m = 0; m < opt.nmodes; ++m)
		{
			//Find mode with highest amplitude
			amplitudes.maxCoeff(&x_mode, &i_mode); //x_mode is always 0 since amplitudes is a 1xN matrix
			BLACS::COMM_ACTIVE.Barrier();
			//cout << BLACS::myrank << " i_mode initial " << i_mode << endl << endl;
			BLACS::COMM_ACTIVE.Barrier();
			// Only print 1 mode out of a possible conjugate pair
			while (lambdas(i_mode, 0).imag() < 0)
			{
				BLACS::COMM_ACTIVE.Barrier();
				amplitudes(0, i_mode) = -1;	// Get the bottom part of the spectrum out of the game
				amplitudes.maxCoeff(&x_mode, &i_mode);
				//cout << BLACS::myrank << " i_mode in loop " << i_mode << endl << endl;
			}
			BLACS::COMM_ACTIVE.Barrier();
			//cout << BLACS::myrank << " i_mode final " << i_mode << endl << endl;
			amplitudes(0, i_mode) = -1;
			BLACS::COMM_ACTIVE.Barrier();
			//cout << BLACS::myrank << " test " << i_mode << endl << endl;

			if (BLACS::mycol == BLACS::indxg2p(i_mode, Modes.cblock(), BLACS::grid_cols))
			{
				int i_loc = BLACS::indxg2l(i_mode, Modes.cblock(), BLACS::grid_cols);
				MatrixXf export;


				// Gather the global column on row 0
				Matrix<MPI::Request, Dynamic, Dynamic> Irecv_requests;
				Matrix<Matrix<float, Dynamic, Dynamic>, Dynamic, 1> RecvBuffer;
				if (BLACS::myrow == 0)
				{
					export.resize(Modes.rows(), 2);
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
					for (int rb = 0; rb < ceil((double)export.rows() / Modes.rblock()); rb++)
					{
						int roffset = Modes.rblock() * floor(rb / BLACS::grid_rows);
						int _nrows = min(Modes.rblock(), export.rows() - rb*Modes.rblock());
						int pr_owner = rb % BLACS::grid_rows;
						export.block(rb*Modes.rblock(), 0, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 0, _nrows, 1);
						export.block(rb*Modes.rblock(), 1, _nrows, 1) = RecvBuffer(pr_owner, 0).block(roffset, 1, _nrows, 1);
					}
				}

				// Print to disk from row 0
				if (BLACS::myrow == 0)
				{
					cout << "(" << BLACS::myrank << ") "<< " writing mode " << m << "...";
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

							values = export.block(offset_gold, 0, dreader.Np, 1);
							gold_print_values(values, georead, pFile);
							
							fclose(pFile);

							stringstream filenameIM;
							filenameIM << "mode" << setfill('0') << setw(6) << m << "." << var << ".ang";

							pFile = fopen((opt.outdir + filenameIM.str()).c_str(), "wb");

							gold_print_header(m, "Angle", var, pFile);

							values = export.block(offset_gold, 1, dreader.Np, 1);
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

		// Write the .case file
		if (BLACS::myrank == 0)
		{
			std::ofstream ofs(opt.outdir + "dmd.case", std::ofstream::out);
			ofs << "FORMAT" << endl
				<< "type: ensight gold" << endl
				<< "GEOMETRY" << endl
				<< "model: dmd.geo" << endl
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
			prof.toc("SaveModes", "\nModes saved in (s): ");
		else
			prof.toc("SaveModes");
		cout << flush;
		BLACS::COMM_ACTIVE.Barrier();
		/////**************************************************************************************************/
		/////*------------------------------      /PRINT SOME MODES TO HDF5      -----------------------------*/
		/////**************************************************************************************************/

		stringstream profile_data;
		profile_data << opt.outdir << argv[0] << "-" << rank << ".yml";

		//if (BLACS::myrank == 0)
		prof.dump(profile_data.str());

		BLACS::COMM_ACTIVE.Barrier();
	} // end IF ACTIVE
	else
	{
		cout << "Pfff" << endl;
	}

	if (ROOT)	// not having if(root) prevents segfault at the begining of the program (!?!)
		cout << endl << endl << "DYMODE OUT!" << endl;

	BLACS::finalize();
	MPI::Finalize();

	return 0;
}

