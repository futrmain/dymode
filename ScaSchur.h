#ifndef PEIGEN_SCASCHUR_H
#define PEIGEN_SCASCHUR_H

#include <boost/math/common_factor.hpp> // For boost::math::lcm


namespace peigen
{
	namespace BLACS
	{		
		extern "C" 
		{
#ifdef _WIN32 /* Win32 or Win64 environment */
#undef numroc_ 
#undef descinit_ 
#undef chk1mat_ 
#undef pchk2mat_ 
#endif 
			int numroc_(int*a, int*b, int*c, int*d, int*e){ return NUMROC(a, b, c, d, e); }
			//void blacs_gridinfo(int* icontxt, int* nprow, int*  npcol, int* myprow, int* mypcol);
			//void blacs_gridinfo_(int* icontxt, int* nprow, int*  npcol, int* myprow, int* mypcol){ cblacs_gridinfo(*icontxt, nprow, npcol, myprow, mypcol); };
			void BLACS_GRIDINFO(int*, int*, int*, int*, int*);
			void blacs_gridinfo_(int*a, int*b, int*c, int*d, int*e) { Cblacs_gridinfo(*a, b, c, d, e); }
			auto descinit_ = DESCINIT;
			void chk1mat_(int*a, int*b, int*c, int*d, int*e, int*f, int*g, int*h, int*i){ CHK1MAT(a, b, c, d, e, f, g, h, i); }
			auto dgsum2d_ = Cdgsum2d;
			void pchk2mat_(int*a, int*b, int*c, int*d, int*e, int*f, int*g, int*h, int*i, int*j, int*k, int*l, int*m, int*n, int*o, int*p, int*q, int*r, int*s, int*t) { PCHK2MAT(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t); }

			void DGEBS2D(MKL_INT *contxt, char *scope, char *top, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *lda);
			auto dgebs2d_ = DGEBS2D;

			//bool LSAME(char *ca, char *cb);
			//auto lsame_ = LSAME;

			void INFOG2L(int* grindx, int* gcindx, int* desc, int* nprow, int* npcol, int*  myrow, int* mycol, int* lrindx, int* lcindx, int* rsrc, int* csrc);
			void infog2l_(int* grindx, int* gcindx, int* desc, int* nprow, int* npcol, int*  myrow, int* mycol, int* lrindx, int* lcindx, int* rsrc, int* csrc)
			{
				INFOG2L( grindx,  gcindx,  desc,  nprow,  npcol,   myrow,  mycol,  lrindx,  lcindx,  rsrc,  csrc);
			}

			void DGEBR2D(int* icontxt, char* scope, char* top, int* m, int* n, double* a, int* lda, int* rsrc, int* csrc);
			auto dgebr2d_ = DGEBR2D;

			void PDELSET(double* a, int* ia, int* ja, int* desca, double* alpha);
			void pdelset_(double* a, int* ia, int* ja, int* desca, double* alpha){ PDELSET(a, ia, ja, desca,alpha); }

			void PDELGET(char* scope, char* top, double* alpha, double* a, int* ia, int* ja, int* desca);
			void pdelget_(char* scope, char* top, double* alpha, double* a, int* ia, int* ja, int* desca){ PDELGET(scope, top, alpha, a, ia, ja, desca); }

			//void dlanv2(double* a, double* b, double* c, double*  d, double* rt1r, double*  rt1i, double* rt2r, double* rt2i, double* cs, double* sn);
			//auto dlanv2_ = dlanv2;

			int ICEIL(int* inum, int* idenom);
			int iceil_(int* inum, int* idenom) { return ICEIL( inum,  idenom);	}

			void DGESD2D(int* icontxt, int* m, int* n, double* a, int* lda, int* rdest, int* cdest);
			void dgesd2d_(int* icontxt, int* m, int* n, double* a, int* lda, int* rdest, int* cdest)
			{ DGESD2D( icontxt,  m,  n,  a,  lda,  rdest,  cdest);	}

			void DGERV2D(int* icontxt, int* m, int*  n, double* a, int* lda, int*  rsrc, int* csrc);
			void dgerv2d_(int* icontxt, int* m, int*  n, double* a, int* lda, int*  rsrc, int* csrc)
			{ DGERV2D( icontxt,  m,   n,  a,  lda,   rsrc,  csrc); }

			void pxerbla(int* ictxt, char srname[6], int* info);
			//auto pxerbla_ = pxerbla;

#ifdef _WIN32 /* Win32 or Win64 environment */
#define numroc_ NUMROC
#define descinit_ DESCINIT
#define chk1mat_ CHK1MAT
#define pchk2mat_ PCHK2MAT
#endif 
		}
	}
}
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

		SharedMatrix<MatrixType>& matrixT() const { return T; }
		SharedMatrix<MatrixType>& matrixZ() const { return Z; }
	};

	template <typename MatrixType>
	ScaSchur<MatrixType>::ScaSchur(const SharedMatrix<MatrixType>& H, const SharedMatrix<MatrixType>& Q, bool computeSchur, bool computeVectors) : T(H), Z(Q), eigenvals(MatrixXcd(H.rows(), 1))
	{
		assert(T.rblock() >= 6 && "Row blocking must be at least 6 rows for PDHSEQR");
		

		char job = computeSchur ? 'S' : 'E';
		char compz = computeVectors ? 'V' : 'N';
		const int N = H.rows();

		double lwork;		
		int info;

		MatrixXd wr(eigenvals.rows(), 1);
		MatrixXd wi(eigenvals.rows(), 1);
	
		int liwork; 

		PBLAS::peigen_pxhseqr(job, compz, N, 1, N, T.localData(), T.desc, wr.data(), wi.data(), Z.localData(), Z.desc, &lwork, -1, &liwork, -1, &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem querying work space for Schur decomposition, return value was " << info << endl;

		//std::cout << "(" << BLACS::myrank << ") SIZE FOR hseqr , info, lwork, liwork: " << info << ", " << lwork << ", " << liwork << std::endl;

		liwork = max((int)lwork, abs(liwork)); // dirty hack because doing a size query sometimes returns negative values for liwork...
		
		cout << BLACS::myrank << " liwork: " << liwork << endl;
		MatrixXd work(liwork, 1);

		
		MatrixXi iwork(liwork, 1);

		cout << "(" << BLACS::myrank << ") " << "Right here, Right now" << endl;
		PBLAS::peigen_pxhseqr(job, compz, N, 1, N, T.localData(), T.desc, wr.data(), wi.data(), Z.localData(), Z.desc, work.data(), (int)lwork, iwork.data(), liwork, &info);
		if (info != 0)
			std::cout << "(" << BLACS::myrank << ") " << "had a problem computing Schur decomposition, return value was " << info << endl;
		else
			std::cout << "(" << BLACS::myrank << ") " << "Schur passed A-OK " << info << endl;

		/*std::cout << "(" << BLACS::myrank << ") SIZE FOR hseqr , info after schur: " << info << std::endl;
		std::cout << "wr " << endl << wr << endl << endl;
		std::cout << "wi " << endl << wi << endl;*/

		for (int k = 0; k < eigenvals.rows(); ++k)
		{
			eigenvals(k, 0) = complex<double>(wr(k), wi(k));
		}
	}

	/**********************************************************************************************************/
	/************************************* OLD CODE USING PDLAHQR *********************************************/
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