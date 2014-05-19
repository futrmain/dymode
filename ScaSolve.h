#ifndef PEIGEN_SCASOLVE_H
#define PEIGEN_SCASOLVE_H

namespace peigen
{
	enum LinSolverName_t { Eigen = 0, pxgesv, pxgesvx };

	template <typename MatrixType>
	class ScaSolve
	{
	public:
		LinSolverName_t method;
		SharedMatrix<MatrixType> solution, matrixLU;

		ScaSolve<MatrixType>(const SharedMatrix<MatrixType>& A, const SharedMatrix<MatrixType>& B, LinSolverName_t solver = pxgesv);

		inline class ScaSolve<MatrixType>& solvername(LinSolverName_t solver)
		{
			method = solver;
			return *this;
		}
	};

	template <typename MatrixType>
	ScaSolve<MatrixType>::ScaSolve(const SharedMatrix<MatrixType>& A, const SharedMatrix<MatrixType>& B, LinSolverName_t solver) : matrixLU(A), solution(B), method(solver)
	{
		assert((A.rblock() == A.cblock()) && "EIGEN SOLVER: 'A' must be distributed using a square block-cyclic distribution");
		assert((A.rblock() == B.rblock()) && "EIGEN SOLVER: 'A' and 'B' must have the same row-block size");

		switch (method)
		{
		case pxgesv:
			int ipiv[matrixLU.local_matrix.rows() + matrixLU.rblock()];	// size described in MKL manual, don't know what +rblock is for
			int desc[9];	// according to IBM doc, desc is modified on exit, and contains info for ipiv instead, so we pass a copy that can be modified safely
			std::copy(matrixLU.descriptor(), matrixLU.descriptor() + 9, desc);

			//A.printDetails();
			//B.printDetails();

			int info;
			PBLAS::pxgesv(matrixLU.x, solution.cols(),
				matrixLU.localData(), matrixLU.i, matrixLU.j, desc,
				ipiv, solution.localData(), solution.i, solution.j, solution.descriptor(), &info);

			if (BLACS::myrank == 0)
				std::cout << "Solved a " << matrixLU.rows() << " x " << matrixLU.cols() << " problem with " << solution.cols() << " rhs. Return code was: " << info << std::endl;
			break;
		}
	}

}	// end namespace peigen
#endif // PEIGEN_SCASOLVE_H
