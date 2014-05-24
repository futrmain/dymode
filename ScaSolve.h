#ifndef PEIGEN_SCASOLVE_H
#define PEIGEN_SCASOLVE_H

namespace peigen
{
	enum LinSolverName_t { Eigen = 0, pxgesv, pxgesvx, svd };

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

		switch (method)
		{
		case pxgesv:
		{
					   assert((A.rblock() == A.cblock()) && "EIGEN SOLVER: 'A' must be distributed using a square block-cyclic distribution");
					   assert((A.rblock() == B.rblock()) && "EIGEN SOLVER: 'A' and 'B' must have the same row-block size");

					   // size described in MKL manual, don't know what +rblock is for
					   int ipiv[matrixLU.local_matrix.rows() + matrixLU.rblock()];
					   ipiv[0] = 0; // initialize to zero to avoid "uninitialized memory access" flooding
					   int desc[9];	// according to IBM doc, desc is modified on exit, and contains info for ipiv instead, so we pass a copy that can be modified safely
					   std::copy(matrixLU.descriptor(), matrixLU.descriptor() + 9, desc);

					   //A.printDetails();
					   //B.printDetails();

					   int info = 0;
					   PBLAS::pxgesv(matrixLU.x, solution.cols(),
						   matrixLU.localData(), matrixLU.i, matrixLU.j, desc,
						   ipiv, solution.localData(), solution.i, solution.j, solution.descriptor(), &info);

					   if (BLACS::myrank == 0)
						   std::cout << "Solved a " << matrixLU.rows() << " x " << matrixLU.cols() << " problem with " << solution.cols() << " rhs. Return code was: " << info << std::endl;

		//A.printDetails();
		//B.printDetails();

					   break;
		}
		case svd:
		{
					matrixLU.gather(0);
					solution.gather(0);

					if (BLACS::ROOT)
					{
						solution.global_matrix = matrixLU.global_matrix.jacobiSvd(ComputeThinU | ComputeThinV).solve(solution.global_matrix);
					}
					solution.dispatch(0);

					break;
		}
		}
	}

}	// end namespace peigen
#endif // PEIGEN_SCASOLVE_H
