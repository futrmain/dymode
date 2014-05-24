#ifndef PEIGEN_SCASCHUR_H
#define PEIGEN_SCASCHUR_H

namespace peigen
{

	template <typename MatrixType>
	class ScaSchur
	{
	protected:
		typedef typename MatrixType::Scalar Scalar;

	public:
		SharedMatrix<MatrixType> T, Z;
		MatrixXcd eigenvals;

		ScaSchur<MatrixType>(SharedMatrix<MatrixType> H, SharedMatrix<MatrixType> Q);
	};

	template <typename MatrixType>
	ScaSchur<MatrixType>::ScaSchur(SharedMatrix<MatrixType> H, SharedMatrix<MatrixType> Q) : T(H), Z(Q)
	{
		assert(T.rows() == T.cols() && "CALLING EIGEN SOLVER ON NON SQUARE MATRIX");


		const int n = T.rows();
		int info;

		int ilo = 1;
		int ihi = n;

		MatrixType work(1, 1);

		MatrixType r(n, 1);
		MatrixType i(n, 1);
		MatrixXi iwork(1, 1);
		iwork(0, 0) = - 42;
		MatrixXi liwork(1, 1);

		int zero = 0;
		int HB = T.rblock();
		int N = n;
		int LOCALK = BLACS::numroc_(&N, &HB, &BLACS::mycol, &zero, &BLACS::grid_cols);
		int JJ = n / T.rblock();


		int LWORK = 3 * n + max(2 * max(*(T.descriptor() + 8), *(T.descriptor() + 8)) + 2 * LOCALK, JJ);
		work.resize(3*LWORK, 1);
		iwork.resize(3*LWORK, 1);

		SharedMatrix<MatrixXd> Z1(Z.rows(), Z.cols(), 'i', Z.rblock(), Z.cblock());

		//Z = Z1;

		//PBLAS::pxlahqr(true, true, n, ilo, ihi, T.localData(), T.descriptor(), r.data(), i.data(), ilo, ihi, Z.localData(), Z.descriptor(), work.data(), LWORK, iwork.data(), LWORK, &info);

		
			int reclevel = 0;

			PBLAS::pxlaqr0(true, true, n, ilo, ihi, T.localData(), T.descriptor(), r.data(), i.data(), ilo, ihi, Z.localData(), Z.descriptor(), work.data(), -1, iwork.data(), -1, &info, reclevel);

			std::cout << "(" << BLACS::myrank << ") SIZE FOR EIGEN SOLVER 1 , info, lwork, iwork: " << info << ", " << work(0, 0) << ", " << iwork(0, 0) << std::endl;


		eigenvals.resize(r.size(), 1);
		eigenvals.real() = r;
		eigenvals.imag() = i;
		
		//std::cout << "(" << BLACS::myrank << ") SIZE FOR EIGEN SOLVER 1 , info, lwork, iwork: " << info << ", " << work(0, 0) << ", " << iwork(0, 0) << std::endl;
		//std::cout << "(" << BLACS::myrank << ") eigenvalues: " << eigenvals << endl << endl << std::endl;

		//if (info != 0)
			std::cout << "(" << BLACS::myrank << ") Schur returned " << info << endl  << std::endl;
	}

}	// end namespace peigen
#endif // PEIGEN_SCASCHUR_H