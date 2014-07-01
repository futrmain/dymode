#ifndef PEIGEN_SCAHESSENBERG_H
#define PEIGEN_SCAHESSENBERG_H

namespace peigen
{

	template <typename MatrixType>
	class ScaHessenberg
	{
	private:
		typedef typename MatrixType::Scalar Scalar;
		void computeQ();
		void computeH();


	public:
		SharedMatrix<MatrixType> Hessenberg;
		SharedMatrix<MatrixType> Q;
		MatrixType tau;
		ScaHessenberg<MatrixType>(const SharedMatrix<MatrixType>& A, const bool computeMatrixQ = true);

		SharedMatrix<MatrixType>&  matrixH() { return Hessenberg; }
		SharedMatrix<MatrixType>&  matrixQ() { return Q; }
	};

	template <typename MatrixType>
	ScaHessenberg<MatrixType>::ScaHessenberg(const SharedMatrix<MatrixType>& A, const bool computeMatrixQ) : Hessenberg(A), tau(MatrixType(A.rows(), 1))
	{
		assert(A.rows() == A.cols() && "CALLING EIGEN SOLVER ON NON SQUARE MATRIX");// FIXME Maybe this is valid for just Hessenberg?

		int info;
		MatrixType work(1, 1);

		PBLAS::pxgehrd(Hessenberg.rows(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(), tau.data(), work.data(), /*lwork*/ -1, &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem querying work space for Hessenberg decomposition, return value was " << info << endl;

		work.resize(work(0, 0), 1);

		PBLAS::pxgehrd(Hessenberg.rows(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(), tau.data(), work.data(), /*lwork*/ work.size(), &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem computing Hessenberg decomposition, return value was " << info << endl; 

		if (computeMatrixQ == true)
		{	// ComputeQ needs to be called before computeH, since the latter destroys information situated in H needed to form Q
			computeQ();
		}

		computeH();
	}

	template <typename MatrixType>
	void ScaHessenberg<MatrixType>::computeH()
	{
		SharedMatrix<MatrixType>& H = Hessenberg;

		//Remove values below subdiagonal, they contain information relative to Q
		for (int r = 2; r < H.rows(); ++r)
		{
			for (int c = 0; c < r - 1; ++c)
			{
				if ((BLACS::myrow == BLACS::indxg2p(r, H.rblock(), BLACS::grid_rows)) && (BLACS::mycol == BLACS::indxg2p(c, H.cblock(), BLACS::grid_cols)))
					H.local_matrix(BLACS::indxg2l(r, H.rblock(), BLACS::grid_rows), BLACS::indxg2l(c, H.cblock(), BLACS::grid_cols)) = 0;
			}
		}
	}

	template <typename MatrixType>
	void ScaHessenberg<MatrixType>::computeQ()
	{
		const int n = Hessenberg.rows();
		int info;

		/*if (BLACS::myrank == 0)
			std::cout << "tau on entry to matrixQ(): " << endl << tau << std::endl << endl;*/

		// Initialize Q with Identity matrix
		Q = SharedMatrix<MatrixType>(n, n, 'i', Hessenberg.rblock(), Hessenberg.cblock());

		MatrixType work(1, 1);

		PBLAS::pxormhr('R', 'N', Q.rows(), Q.cols(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(),
			tau.data(), Q.localData(), Q.i, Q.j, Q.descriptor(), work.data(), -1, &info);

		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem querying work space to compute Q from the Hessenberg decomposition, return value was " << info << endl;

		//if (BLACS::myrank == 0)
		//	std::cout << "SIZE FOR ORMHR, lwork: " << work(0, 0) << std::endl;
		work.resize(work(0, 0), 1);

		PBLAS::pxormhr('R', 'N', Q.rows(), Q.cols(), /*ilo*/ Hessenberg.x, /*ihi*/ Hessenberg.y, Hessenberg.localData(), Hessenberg.i, Hessenberg.j, Hessenberg.descriptor(),
			tau.data(), Q.localData(), Q.i, Q.j, Q.descriptor(), work.data(), work.size(), &info);
		
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem computing Q from the Hessenberg decomposition, return value was " << info << endl;
	}


}	// end namespace peigen
#endif // PEIGEN_SCAHESSENBERG_H
