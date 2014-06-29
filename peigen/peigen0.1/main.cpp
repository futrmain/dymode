#include <iostream>
#include <Eigen/Dense>
#include <fstream>

#include "mpi.h"

#include "BLACS.h"
#include "PBLAS.h"
#include "SharedMatrix.h"

using namespace std;
using namespace Eigen;
using namespace peigen;

int main(int argc, char* argv[])
{
	MPI::Init();

	int rank, numtasks;
	rank = MPI::COMM_WORLD.Get_rank();
	numtasks = MPI::COMM_WORLD.Get_size();
	bool ROOT = (rank==0);

	// create the BLACS grid
	BLACS::init(numtasks);


	SharedMatrix<MatrixXd> A = SharedMatrix<MatrixXd>::Proc(10,8);
	SharedMatrix<MatrixXd> B = SharedMatrix<MatrixXd>::Proc(6,5, 5,3);

	if(ROOT)
		cout << "A: " << endl;
	cout << A;
	BLACS::Cblacs_barrier(BLACS::ctxt, "All");

	if(ROOT)
		cout << "B: " << endl;
	cout << B;
	BLACS::Cblacs_barrier(BLACS::ctxt, "All");

	//SharedMatrix<MatrixXd> P = A.block(1,2, 5,6).transpose() * B.transpose().block(0,0, 5,4);

	/*if(ROOT)
		cout << "P: " << P.rows() << "x"<<P.cols()<< flush <<  endl;
	cout << P;
	BLACS::Cblacs_barrier(BLACS::ctxt, "All");*/

	
	SharedMatrix<MatrixXd> C = SharedMatrix<MatrixXd>::Ones(4, 9);

	SharedMatrix<MatrixXd> P2 = A.block(1,2, 5,6).transpose() * B.transpose().block(0,0, 5,4) * C;

	if(ROOT)
		cout << "P2: " << P2.rows() << "x"<<P2.cols()<< flush <<  endl;
	cout << P2;
	BLACS::Cblacs_barrier(BLACS::ctxt, "All");

	if(ROOT)
		cout << "B: " << B.rows() << "x"<<B.cols()<< flush <<  endl;
	if(ROOT)
		cout << "B: " << B.x << "x"<<B.y<< flush <<  endl;
	cout << B;
	BLACS::Cblacs_barrier(BLACS::ctxt, "All");

	MPI::Finalize();		    
	return 0;
}