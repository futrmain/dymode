#ifndef PEIGEN_SCAHESSENBERG_H
#define PEIGEN_SCAHESSENBERG_H

namespace peigen
{

	template <typename MatrixType>
	class ScaHessenberg
	{
	protected:
		typedef typename MatrixType::Scalar Scalar;

	public:
		SharedMatrix<MatrixType> H, Q;

		ScaHessenberg<MatrixType>(SharedMatrix<MatrixType> A);
	};

	template <typename MatrixType>
	ScaHessenberg<MatrixType>::ScaHessenberg(SharedMatrix<MatrixType> A) : Q(SharedMatrix<MatrixType>(A.rows(), A.rows(), 'i', A.rblock(), A.cblock()))
	{
		assert(A.rows() == A.cols() && "CALLING EIGEN SOLVER ON NON SQUARE MATRIX");

		/*Step 1: Form H*/

		H = A;
		const int n = H.rows();
		int info;

		int ilo = 1;
		int ihi = n;

		MatrixType work(1, 1);
		MatrixType tau(n, 1);

		PBLAS::pxgehrd(n, /*ilo*/ ilo, /*ihi*/ ihi, H.localData(), H.i, H.j, H.descriptor(), tau.data(), work.data(), /*lwork*/ -1, &info);

		//if (BLACS::myrank == 0)
		//std::cout << "SIZE FOR HESSENBERG, lwork: " << work(0, 0) << std::endl;

		work.resize(work(0, 0), 1);
		//tau.resize(tau(0, 0), 1);

		PBLAS::pxgehrd(n, /*ilo*/ ilo, /*ihi*/ ihi, H.localData(), H.i, H.j, H.descriptor(), tau.data(), work.data(), /*lwork*/ work.size(), &info);

		if (info!=0)
			std::cout << "(" << BLACS::myrank << ") Hessenberg returned " << info << std::endl;



		/*Step 2: Form Q*/

		//Q = SharedMatrix<MatrixType>(n, n, 'i', H.rblock(), H.cblock());
		PBLAS::pxormhr('R', 'N', n, n, 1, n, H.localData(), H.i, H.j, H.descriptor(), tau.data(), Q.localData(), Q.i, Q.j, Q.descriptor(), work.data(), -1, &info);

		if (BLACS::myrank == 0)
			std::cout << "SIZE FOR ORMHR, lwork: " << work(0, 0) << std::endl;
		work.resize(work(0, 0), 1);

		PBLAS::pxormhr('R', 'N', n, n, 1, n, H.localData(), H.i, H.j, H.descriptor(), tau.data(), Q.localData(), Q.i, Q.j, Q.descriptor(), work.data(), work.size(), &info);


		if (info != 0)
		std::cout << "(" << BLACS::myrank << ") Q returned " << info << std::endl;

		/*Step 3: Remove values below subdiagonal, they have been used to form Q*/
		for (int r = 2; r < H.rows(); ++r)
			for (int c = 0; c < r - 1; ++c)
				if ((BLACS::myrow == BLACS::indxg2p(r, H.rblock(), BLACS::grid_rows)) && (BLACS::mycol == BLACS::indxg2p(c, H.cblock(), BLACS::grid_cols)))
					H.local_matrix(BLACS::indxg2l(r, H.rblock(), BLACS::grid_rows), BLACS::indxg2l(c, H.cblock(), BLACS::grid_cols)) = 0;
	}

}	// end namespace peigen
#endif // PEIGEN_SCAHESSENBERG_H