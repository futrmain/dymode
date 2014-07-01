#ifndef PEIGEN_SCASOLVE_H
#define PEIGEN_SCASOLVE_H

#include <type_traits>

namespace peigen
{
	enum LinSolverName_t { Eigen = 0, pxgesv, pxgesvx, EigenSVD };

	template< typename MatrixType, bool isComplex = (is_same<typename MatrixType::Scalar, complex<float>>::value || is_same<typename MatrixType::Scalar, complex<double>>::value) >
	class xwork 
	{
	private:
		Matrix<int, Dynamic, 1> work;
	public:
		int* data() { return work.data(); }
		int operator()(int i, int j) { return work(i, j); }
		void resize(int x, int y) { work.resize(x, y); }
		int rows() { return work.rows(); }
		xwork() : work(Matrix<int, Dynamic, 1>(1, 1)) {}

	};

	template< typename MatrixType >
	class xwork<MatrixType, true >
	{
	private:
		typedef typename MatrixType::RealScalar RealScalar;
		typedef typename MatrixType::Scalar Scalar;
		Matrix<RealScalar, Dynamic, 1> work;
	public:
		Scalar* data() { return work.data(); }
		Scalar operator()(int i, int j) { return work(i, j); }
		void resize(int x, int y) { work.resize(x, y); }

		xwork() : work(Matrix<Scalar, Dynamic, 1>(1, 1)) {}
	};

	template <typename MatrixType>
	class ScaSolve
	{
	private:
	  typedef typename MatrixType::RealScalar RealScalar;
	public:
		LinSolverName_t method;
		SharedMatrix<MatrixType> solution, matrixLU;

		ScaSolve<MatrixType>(const SharedMatrix<MatrixType>& A, const SharedMatrix<MatrixType>& B, LinSolverName_t solver = pxgesv);

		inline class ScaSolve<MatrixType>& solvername(LinSolverName_t solver)
		{
			method = solver;
			return *this;
		}

		RealScalar residual(SharedMatrix<MatrixType> A, SharedMatrix<MatrixType> B)
		{
			SharedMatrix<MatrixType> R(B);
			R.pgemm(1., A, solution, -1.);
			
			// Compute local highest residual
			return R.local_matrix.cols() > 0 ? R.local_matrix.cwiseAbs().maxCoeff() : -1;
		}

		RealScalar global_residual(SharedMatrix<MatrixType> original, SharedMatrix<MatrixType> rhs)
		{
			double r_loc = residual(original, rhs);
			double r;

			BLACS::COMM_ACTIVE.Allreduce(&r_loc, &r, 1, MPI::DOUBLE, MPI::MAX);

			return r;
		}
	};

	template <typename MatrixType>
	ScaSolve<MatrixType>::ScaSolve(const SharedMatrix<MatrixType>& A, const SharedMatrix<MatrixType>& B, LinSolverName_t solver) : matrixLU(A), solution(B), method(solver)
	{

		switch (method)
		{
		// FIXME missing implementations for comparison
		case Eigen:
		{
					  assert(0 && "LINEAR SOLVER: Solving with Eigen has nor been implemented yet.");
					  break;
		}
		case pxgesvx:
		{
						assert((A.rblock() == A.cblock()) && "LINEAR SOLVER: 'A' must be distributed using a square block-cyclic distribution");
						assert((A.rblock() == B.rblock()) && "LINEAR SOLVER: 'A' and 'B' must have the same row-block size");

						// size described in MKL manual, don't know what +rblock is for
						int ipiv[matrixLU.local_matrix.rows() + matrixLU.rblock()];
						ipiv[0] = 0; // initialize to zero to avoid "uninitialized memory access" flooding
						int desc[9];	// according to IBM doc, desc is modified on exit, and contains info for ipiv instead, so we pass a copy that can be modified safely
						std::copy(matrixLU.descriptor(), matrixLU.descriptor() + 9, desc);

						//A.printDetails();
						//B.printDetails();

						int info = 0;

						SharedMatrix<MatrixType> Af(A);
						SharedMatrix<MatrixType> x(B);
						Matrix<RealScalar, Dynamic, 1> r(matrixLU.local_matrix.rows(), 1);
						Matrix<RealScalar, Dynamic, 1> c(matrixLU.local_matrix.cols(), 1);

						RealScalar rcond;
						Matrix<RealScalar, Dynamic, 1> ferr(solution.local_matrix.cols(), 1);
						Matrix<RealScalar, Dynamic, 1> berr(solution.local_matrix.cols(), 1);

						MatrixType work(1, 1);
						xwork<MatrixType> xspace;


						PBLAS::pxgesvx('E', 'N', matrixLU.x, solution.cols(),
							matrixLU.localData(), matrixLU.i, matrixLU.j, matrixLU.desc,
							Af.localData(), Af.i, Af.j, Af.desc,
							ipiv, 'N', r.data(), c.data(),
							solution.localData(), solution.i, solution.j, solution.desc,
							x.localData(), x.i, x.j, x.desc,
							&rcond, ferr.data(), berr.data(),
							work.data(), -1, xspace.data(), -1, &info);

						work.resize((int)work(0, 0), 1);
						xspace.resize((int)xspace(0, 0), 1);

						PBLAS::pxgesvx('E', 'N', matrixLU.x, solution.cols(),
							matrixLU.localData(), matrixLU.i, matrixLU.j, matrixLU.desc,
							Af.localData(), Af.i, Af.j, Af.desc,
							ipiv, 'N', r.data(), c.data(),
							solution.localData(), solution.i, solution.j, solution.desc,
							x.localData(), x.i, x.j, x.desc,
							&rcond, ferr.data(), berr.data(),
							work.data(), work.rows(), xspace.data(), xspace.rows(), &info);


						/*if (BLACS::myrank == 0)
							std::cout << "lwork is " << work(0, 0) << ", lrwork is: " << xspace(0,0) << ", info is: " << info << std::endl;*/

						if (BLACS::myrank == 0)
							std::cout << "rcond = " << rcond << endl << "ferr = " << ferr(0, 0) << endl << "berr = " << berr(0,0) << std::endl;

						if (BLACS::myrank == 0)
							std::cout << "Solved a " << matrixLU.rows() << " x " << matrixLU.cols() << " problem with " << solution.cols() << " rhs using pxgesvx." << endl << "Return code was: " << info << endl << std::endl;

						solution = x;

						break;
		}
		case pxgesv:
		{
					   assert((A.rblock() == A.cblock()) && "LINEAR SOLVER: 'A' must be distributed using a square block-cyclic distribution");
					   assert((A.rblock() == B.rblock()) && "LINEAR SOLVER: 'A' and 'B' must have the same row-block size");

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
						   std::cout << "Solved a " << matrixLU.rows() << " x " << matrixLU.cols() << " problem with " << solution.cols() << " rhs using pxgesvx." << endl << "Return code was: " << info << endl << std::endl;

		//A.printDetails();
		//B.printDetails();

					   break;
		}
		case EigenSVD:
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
