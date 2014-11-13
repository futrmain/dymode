#ifndef PEIGEN_SCASCHUR_H
#define PEIGEN_SCASCHUR_H

#include <boost/math/common_factor.hpp> // For boost::math::lcm

namespace peigen
{

	template <typename MatrixType>
	class ScaSchur
	{
	protected:
		typedef typename MatrixType::Scalar Scalar;

	private:
		SharedMatrix<MatrixType> T, Z;

	public:
		MatrixXcd eigenvals;

		ScaSchur<MatrixType>(const SharedMatrix<MatrixType>& H, const SharedMatrix<MatrixType>& Q, bool computeSchur = true, bool compz = true);

		SharedMatrix<MatrixType>& matrixT() { return T; }
		SharedMatrix<MatrixType>& matrixZ() { return Z; }
	};

	template <typename MatrixType>
	ScaSchur<MatrixType>::ScaSchur(const SharedMatrix<MatrixType>& H, const SharedMatrix<MatrixType>& Q, bool computeSchur, bool computeVectors) : T(H), Z(Q), eigenvals(MatrixXcd(H.rows(), 1))
	{
		if (ROOT)
		{
			std::cout << "Warning: the ScaLAPACK routine p*hseqr used for Schur reduction may not work in all ScaLAPACK implementation." << std::endl;
		}
		assert(T.rblock() >= 6 && "Row blocking must be at least 6 rows for PDHSEQR");
		

		char job = computeSchur ? 'S' : 'E';
		char compz = computeVectors ? 'V' : 'N';
		const int N = H.rows();

		double lwork;		
		int info;

		MatrixXd wr(eigenvals.rows(), 1);
		MatrixXd wi(eigenvals.rows(), 1);
	
		int liwork; 

		PBLAS::pxhseqr(job, compz, N, 1, N, T.localData(), T.desc, wr.data(), wi.data(), Z.localData(), Z.desc, &lwork, -1, &liwork, -1, &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem querying work space for Schur decomposition, return value was " << info << endl;

		//std::cout << "(" << BLACS::myrank << ") SIZE FOR hseqr , info, lwork, liwork: " << info << ", " << lwork << ", " << liwork << std::endl;

		liwork = max((int)lwork, abs(liwork)); // dirty hack because doing a size query sometimes returns negative values for liwork...
		// On 2014-10-10 MKL dev team said the next update of MKL will have fixed the problems with Schur decomposition.

		liwork = max(1, liwork);	// As recomended on Intel forum. Probably a best practice since ScaLAPACK sometimes returns 0 as a work size but still accesses the adress given. 
		
		//cout << BLACS::myrank << " liwork: " << liwork << endl;
		MatrixXd work(liwork, 1);

		
		MatrixXi iwork(liwork, 1);

		//cout << "(" << BLACS::myrank << ") " << "Right here, Right now" << endl;
		PBLAS::peigen_pxhseqr(job, compz, N, 1, N, T.localData(), T.desc, wr.data(), wi.data(), Z.localData(), Z.desc, work.data(), (int)lwork, iwork.data(), liwork, &info);
		PBLAS::pxhseqr(job, compz, N, 1, N, T.localData(), T.desc, wr.data(), wi.data(), Z.localData(), Z.desc, work.data(), (int)lwork, iwork.data(), liwork, &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem computing Schur decomposition, return value was " << info << endl;
		else
			std::cout << "(" << BLACS::myrank << ") " << "Schur passed A-OK, returned info was " << info << endl;

		/*std::cout << "(" << BLACS::myrank << ") SIZE FOR hseqr , info after schur: " << info << std::endl;
		std::cout << "wr " << endl << wr << endl << endl;
		std::cout << "wi " << endl << wi << endl;*/

		for (int k = 0; k < eigenvals.rows(); ++k)
		{
			eigenvals(k, 0) = complex<double>(wr(k), wi(k));
		}
	}

	/**********************************************************************************************************/
	/************************************* OLD CODE USING PDLAHQR FOLLOWS *************************************/
	/**********************************************************************************************************/


	// Kept as comment in case it needs to be implemented sometime
	// PDLAHQR gives a funny Schur matrix, didn't investigate much
	// Also, using -1 to query work size does not work, unlike what MKL's doc says


	//assert(T.rblock() >= 5 && "Row blocking must be at least 5 rows for PDLAHQR");

	//bool options = true;
	//char job = computeSchur ? 'S' : 'E';
	//char compz = 'V'; // Left multiply Z by Q
	//const int N = H.rows();

	//int lwork;		
	//int info;

	//MatrixXd wr(eigenvals.rows(), 1);
	//MatrixXd wi(eigenvals.rows(), 1);

	//lwork = 3 * T.rows() + max(2 * max(Z.desc[9 - 1], T.desc[9 - 1]) + 2 * max(T.local_matrix.rows(), T.local_matrix.cols()), 7 * ceil(T.rows() / T.rblock()) / boost::math::lcm(BLACS::grid_rows, BLACS::grid_cols));

	//int iwork; // iwork is just a placeholder in this version of ScaLAPACK
	//MatrixXd work(lwork, 1);

	//PBLAS::pxlahqr(options, options, T.rows(), 1, T.rows(), T.localData(), T.desc, wr.data(), wi.data(), 1, Z.rows(), Z.localData(), Z.desc, work.data(), lwork, &iwork, -1, &info);
	//
	//std::cout << "(" << BLACS::myrank << ") SIZE FOR EIGEN SOLVER 1 , info, lwork, iwork: " << info << ", " << lwork << ", " << iwork << std::endl;

	//std::cout << "wr " << endl << wr << endl << endl;
	//std::cout << "wi " << endl << wi << endl;

}	// end namespace peigen
#endif // PEIGEN_SCASCHUR_H
