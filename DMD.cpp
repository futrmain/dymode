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


		/*assert(svd.singularValues(0, 0) == 5564.1186373992377 && "Singular value 0 changed");
		assert(svd.singularValues(1, 0) == 22.886435861737159 && "Singular value 1 changed");
		assert(svd.singularValues(2, 0) == 1.8529983756209095e-012 && "Singular value 2 changed");
		assert(svd.singularValues(3, 0) == 8.8044075548124019e-013 && "Singular value 3 changed");
		assert(svd.singularValues(4, 0) == 2.3554640784338274e-013 && "Singular value 4 changed");
		assert(svd.singularValues(5, 0) == 1.2286640032115696e-013 && "Singular value 5 changed");*/

		if (ROOT)
			cout << "HERE COMES U" << endl << flush;
		cout << svd.matrixU;
		svd.matrixU.gather(0);
		if (ROOT)
		{
			MatrixXd Mtest(12, 2);
			Mtest << -0.0033166325649726103, -0.54301086068349469,
				-0.04733862156291422, -0.47191183482749016,
				-0.091360610560856204, -0.40081280897148758,
				-0.13538259955879819, -0.32971378311548916,
				-0.17940458855674019, -0.25861475725948946,
				-0.22342657755468212, -0.18751573140349323,
				-0.26744856655262417, -0.11641670554749042,
				-0.31147055555056613, -0.045317679691494445,
				-0.35549254454850815, 0.02578134616450567,
				-0.39951453354645011, 0.09688037202050867,
				-0.44353652254439202, 0.16797939787650237,
				-0.48755851154233409, 0.23907842373250082;
			/*Mtest <<
				-0.00331663256497242, 0.543010860683501,
				-0.0473386215629142, 0.471911834827486,
				-0.0913606105608562, 0.400812808971487,
				-0.135382599558798, 0.329713783115488,
				-0.17940458855674, 0.258614757259489,
				-0.223426577554682, 0.187515731403491,
				-0.267448566552624, 0.116416705547492,
				-0.311470555550566, 0.045317679691494,
				-0.355492544548508, -0.0257813461645077,
				-0.39951453354645, -0.0968803720205034,
				-0.443536522544392, -0.167979397876503,
				-0.487558511542334, -0.239078423732501;
				Mtest.col(1).noalias() = - Mtest.col(1);*/
			//assert(svd.matrixU.global_matrix.block(0, 0, 12, 2) == Mtest && "Matrix U is different");
		}
		svd.matrixU.global_matrix.resize(0, 0);

		if (ROOT)
			cout << "HERE COMES Vt" << endl << flush;
		cout << svd.matrixVt;
		svd.matrixVt.gather(0);
		if (ROOT)
		{
			MatrixXd Mtest(2, 6);
			Mtest << -0.40426931286211548, -0.40585730067649162, -0.40744528849086736, -0.40903327630524311, -0.41062126411961897, -0.41220925193399488,
				0.60031312370089751, 0.36127267641724625, 0.12223222913361902, -0.11680821815001126, -0.35584866543362653, -0.59488911271727174;
			//assert(svd.matrixVt.global_matrix.block(0, 0, 2, 6) == Mtest && "Matrix Vt is different");
		}
		svd.matrixVt.global_matrix.resize(0, 0);

		svd.matrixU.gather(0);
		svd.matrixVt.gather(0);
		snaps.gather(0);
		if (ROOT)
		{
			MatrixXd Mtest = svd.matrixU.global_matrix * svd.singularValues.asDiagonal() * svd.matrixVt.global_matrix - snaps.global_matrix.block(0, 0, 12, 6);
			cout << "The SVD residuals are " << endl << flush;
			cout << Mtest << endl << endl;
		}
		svd.matrixVt.global_matrix.resize(0, 0);
		svd.matrixU.global_matrix.resize(0, 0);
		snaps.global_matrix.resize(0, 0);

		prof.toc("SVD");
		/////**************************************************************************************************/
		/////*-------------------------------------      /DO AN SVD      -------------------------------------*/
		/////**************************************************************************************************/

		BLACS::COMM_ACTIVE.Barrier();


		///////**************************************************************************************************/
		///////*------------------------------      DO B := Ut * S2 * V * SIG+     ------------------------------*/
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

		B.gather(0);
		if (ROOT)
		{
			MatrixXd Mtest(6, 6);
			Mtest <<
				5585.7612451895984, -0.14377565735435383, 1.4210854715202004e-013, 5.1159076974727213e-013, -2.2737367544323206e-013, -4.5474735088646412e-013,
				13.400333275502359, 22.797415086416748, 1.8829382497642655e-013, 2.9443114613059151e-013, -1.829647544582258e-013, 2.5401902803423582e-013,
				4.7910308660233113e-013, 6.9855291876648023e-014, -2.8414910809283506e-015, -2.6914589092914937e-014, 3.3521899611819412e-014, -5.3986622902111564e-015,
				4.8786164695942578e-013, 7.8289904671112035e-014, 4.2797632188579284e-015, -2.1093230769978214e-015, -3.605445225763021e-014, -1.0404780615398126e-014,
				-2.2096324278813057e-013, -1.8431204655401709e-013, -1.2558665611118508e-027, -2.5306657839490059e-027, 1.7796702021786026e-027, -1.4136387421560532e-027,
				6.953590958353444e-013, -1.2691567635027572e-013, -2.0296004253972114e-014, -5.4354109428787505e-014, -1.5238995944984674e-014, 2.8858286716652247e-014;

			//assert(B.global_matrix == Mtest && "Matrix B partial is different");
		}
		B.global_matrix.resize(0, 0);


		//// Now only a right-multiply by SIG+ is left to do
		//// Each process owns the entire vector sig. Figuring out which value goes to which process
		//// and creating a shared matrix and multiplying is a pain. 
		//// Since SIG+ is diagonal we do what it does: scale the columns of B by the propre value.
		MatrixXd SIGplus = MatrixXd::Zero(B.local_matrix.cols(), B.local_matrix.cols());
		double pinv_tol = std::numeric_limits<double>::epsilon() * (4 * Np) * svd.singularValues(0);
		if (ROOT)
		{
			cout << " pinv_tol: " << pinv_tol << endl << flush;
			assert(pinv_tol == 1.4825790295167806e-011 && "Tolerance for pseudo-inverse differs");
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
			MatrixXd Mtest(1, 6);
			Mtest << 0.00017972298312952879, 0.043694003122253591, 0, 0, 0, 0;
			//assert(SIGplus.diagonal().transpose() == Mtest && "Tolerance for pseudo-inverse differs");
		}


		B.local_matrix = B.local_matrix * SIGplus;
		if (ROOT)
			cout << "B is computed " << endl;
		if (ROOT)
			cout << "HERE COMES Ut M V SIG+" << endl << flush;
		cout << B;

		B.gather(0);
		if (ROOT)
		{
			MatrixXd Mtest(6, 6);
			Mtest <<
				5585.7612451896011, -0.1437756573545812, -1.3642420526593924e-012, 1.9895196601282805e-013, -1.1368683772161603e-012, 6.2527760746888816e-013,
				13.400333275503115, 22.797415086416777, 2.9665159217984183e-013, 3.0531133177191805e-013, 5.3290705182007514e-015, -1.3300471835009375e-013,
				5.7512003390063787e-013, 1.1337998113839656e-013, 3.5691895184709159e-014, -3.513115980210724e-014, -1.0862217136093942e-013, -3.0615889943037428e-014,
				7.1928694561196024e-013, -1.8366825857451476e-014, -1.0437288326033194e-013, -5.0436929207818793e-014, -5.8378645735204742e-014, -1.4235239117574029e-014,
				2.8977808992629942e-013, -4.6087593677805591e-014, 7.2374793893705303e-015, -3.5808371967764976e-014, -6.4648175842202785e-015, -3.7221904766671816e-014,
				0, 0, 0, 0, 0, 0;

			//assert(B.global_matrix == Mtest && "Matrix B partial is different");
		}
		B.global_matrix.resize(0, 0);

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

		// Matrix of eigen vectors
		SharedMatrix<MatrixXcd> X;
		MatrixXcd lambdas(Nt - 1, 1);

		ScaEigenSolver<MatrixXd> eig(B, true, EigSchur);

		double r = eig.global_residual(B);
		if (ROOT)
			cout << "Residual from Eigen problem: " << r << endl << flush;

		prof.toc("EigenProblem");
		/////**************************************************************************************************/
		/////*------------------------------    /DO AN EIGENVECTOR SOLUTION     ------------------------------*/
		/////**************************************************************************************************/

		BLACS::COMM_ACTIVE.Barrier();


		/////**************************************************************************************************/
		/////*------------------------------      DO A LINEAR SYSTEM SOLVE      ------------------------------*/
		/////**************************************************************************************************/
		//prof.tic("LinearSolve");
		////cout << "(" << BLACS::myrank << ")" << endl;
		//SharedMatrix<MatrixXd> rhs = svd.matrixU.transpose() * snaps.block(0, Nt - 1, 4 * Np, 1);
		//svd.matrixU.clear();
		//snaps.clear();
		//SharedMatrix<MatrixXcd> rhsZ = rhs.cast<std::complex<double> >();

		//BLACS::COMM_ACTIVE.Barrier();

		//if (ROOT)
		//	cout << "rhs is computed" << endl << flush;
		//cout << rhsZ;

		//rhsZ.gather(0);
		//if (ROOT)
		//{
		//	MatrixXcd Mtest(6, 1);
		//	Mtest <<
		//		complex<double>(-2302.4169337882709, 0),
		//		complex<double>(-19.085675388167118, 0),
		//		complex<double>(-2.2737367544323206e-013, 0),
		//		complex<double>(-2.8421709430404007e-013, 0),
		//		complex<double>(-8.5265128291212022e-014, 0),
		//		complex<double>(-5.6843418860808015e-014, 0);

		//	assert(rhsZ.global_matrix == Mtest && "Matrix RHS is different");
		//}
		//rhsZ.global_matrix.resize(0, 0);

		//ScaSolve<MatrixXcd> solver(X, rhsZ, peigen::svd);
		//if (ROOT)
		//{
		//	cout << "HERE COMES solution is computed" << endl << flush;
		//}
		//cout << solver.solution << endl;

		//solver.solution.gather(0);
		//if (ROOT)
		//{
		//	MatrixXcd Mtest(6, 1);
		//	Mtest <<
		//		complex<double>(1096070658.1017988, -1322.9523436264617),
		//		complex<double>(1096070658.1017988, 1322.9523435911988),
		//		complex<double>(1.3694775963610608e-011, -7.9879698252414504e-029),
		//		complex<double>(-9.2969462508485929e-012, 5.540911226960229e-029),
		//		complex<double>(7.2585917425865843e-012, -4.3740644386348342e-029),
		//		complex<double>(-5.6843418860808015e-014, 0);

		//	assert(solver.solution.global_matrix == Mtest && "Weights solution is different");
		//}
		//solver.solution.global_matrix.resize(0, 0);
		//

		//if (ROOT)
		//{
		//	cout << "X * w" << endl << flush;
		//}
		//SharedMatrix<MatrixXcd> p1 = X * solver.solution;
		//cout << p1 << endl;

		//if (ROOT)
		//{
		//	cout << "U' * Slast" << endl << flush;
		//}
		//SharedMatrix<MatrixXd> p2 = svd.matrixU.transpose() * snaps.block(0, Nt - 1, 4 * Np, 1);
		//svd.matrixU.clear();
		//snaps.clear();
		//cout << p2 << endl;


		//if (ROOT)
		//{
		//	cout << "RESIDUALS" << endl << flush;
		//}
		//SharedMatrix<MatrixXcd> r = p1;
		//r.local_matrix -= p2.cast<std::complex<double> >().local_matrix;
		//cout << r << endl;
		//


		//SharedMatrix<MatrixXcd> Modes = svd.matrixU.cast<std::complex<double> >() * X;

		//

		//// Each process gets the full, global weighting vector
		//for (int proc = 0; proc < BLACS::grid_rows*BLACS::grid_cols; proc++)
		//{
		//	solver.solution.gather(proc);
		//}

		///*if (ROOT)
		//	cout << "Details for solution " << solver.solution.global_matrix.rows() << endl << flush;
		//	*/
		//BLACS::COMM_ACTIVE.Barrier();

		///*
		//if (ROOT)
		//	cout << "Details for Modes " << endl << flush;
		//Modes.printDetails();
		//*/
		//BLACS::COMM_ACTIVE.Barrier();


		//MatrixXcd weight = MatrixXcd::Zero(Modes.local_matrix.cols(), Modes.local_matrix.cols());
		//for (int i = 0; i < Modes.local_matrix.cols(); i++)
		//{
		//	int index = i % Modes.cblock() + (floor(i / Modes.cblock())*BLACS::grid_cols + BLACS::mycol) * Modes.cblock();
		//	//cout << "(" << BLACS::myrank << ")" << i << " " << index << endl;
		//	weight(i, i) = solver.solution.global_matrix(index);
		//}

		////cout << "(" << BLACS::myrank << ")" << endl;

		//Modes.local_matrix = Modes.local_matrix * weight;

		//if (ROOT)
		//{
		//	cout << "HERE COME the modes" << endl << flush;
		//}
		//cout << Modes;

		//Modes.clear();

		//if (ROOT)
		//	cout << "HERE COMES Modes... Nah just kidding, Modes is just way too big" << endl << flush;
		////cout << Modes;
		//prof.toc("LinearSolve");
		/////**************************************************************************************************/
		/////*------------------------------      /DO A LINEAR SYSTEM SOLVE      -----------------------------*/
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

