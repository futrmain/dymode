#ifndef PEIGEN_SCASOLVE_H
#define PEIGEN_SCASOLVE_H

namespace peigen
{
	template <typename MatrixType>
	class ScaSolve
	{
	public:
		SharedMatrix<MatrixType> solution, matrixLU;

		ScaSolve<MatrixType>(SharedMatrix<MatrixType> A, SharedMatrix<MatrixType> B);
	};

	template <typename MatrixType>
	ScaSolve<MatrixType>::ScaSolve<MatrixType>(SharedMatrix<MatrixType> A, SharedMatrix<MatrixType> B) : matrixLU(A), solution(B) 
	{
		assert((A.rblock() == A.cblock()) && "EIGEN SOLVER: 'A' must be distributed using a square block-cyclic distribution");
		assert((A.rblock() == B.rblock()) && "EIGEN SOLVER: 'A' and 'B' must have the same row-block size");

		int ipiv[matrixLU.local_matrix.rows() + matrixLU.rblock()];	// size described in MKL manual, don't know what +rblock is for
		int desc[9];	// according to IBM doc, desc is modified on exit, and contains info for ipiv instead, so we pass a copy that can be modified safely
		

		std::copy(matrixLU.descriptor(), matrixLU.descriptor()+9, desc);

		A.printDetails();
		B.printDetails();

		int info;
		PBLAS::pxgesv(matrixLU.x, solution.cols(), 
			matrixLU.localData(), matrixLU.i, matrixLU.j, desc,
			ipiv, solution.localData(), solution.i, solution.j, solution.descriptor(), &info);

		if(BLACS::myrank==0)
			std::cout << "Solved a " << matrixLU.rows() << " x " << matrixLU.cols() << " problem with " << solution.cols() << " rhs. Return code was: " << info << std::endl;
	}

}	// end namespace peigen
#endif // PEIGEN_SCASOLVE_H
