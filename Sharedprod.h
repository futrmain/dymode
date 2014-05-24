#ifndef PEIGEN_SHARED_PROD_H
#define PEIGEN_SHARED_PROD_H

namespace peigen
{
	using namespace Eigen;

	template <typename MatrixType>
	class SharedMatrix;

	template <typename MatrixType>
	class Sharedprod
	{
	public:
		typedef typename MatrixType::Scalar Scalar;
		Scalar alpha;
		SharedMatrix<MatrixType> &A;
		SharedMatrix<MatrixType> &B;
		
		Sharedprod<MatrixType>(SharedMatrix<MatrixType> &a, SharedMatrix<MatrixType> &b) : A(a), B(b), alpha(1) {}
		Sharedprod<MatrixType>(Scalar x, Sharedprod<MatrixType> &P) : A(P.A), B(P.B), alpha(x*P.alpha) {}
		Sharedprod<MatrixType>(Sharedprod<MatrixType> &P, SharedMatrix<MatrixType> &b) : A(P.eval()), B(b), alpha(1) {std::cout << A;}

		SharedMatrix<MatrixType> eval();
	};

	template <typename MatrixType>
	SharedMatrix<MatrixType> Sharedprod<MatrixType>::eval()
	{
		// Allocate result matrix
		int m, n, k, rb, cb;

		char opA = (A.use_transpose? 'T' : 'N');	// FIXME need to distinguish transpose and adjoint
		char opB = (B.use_transpose? 'T' : 'N');

		rb = A.rblock();
		cb = B.cblock();

		SharedMatrix<MatrixType> P(A.x, B.y, rb, cb);

		//A.printDetails();

		//B.printDetails();

		//P.printDetails();

		assert((A.y==B.x) && "Multiplying matrices of different size");
		

		// Define submatrices
		int iA = A.i;
		int jA = A.j;
		int iB = B.i;
		int jB = B.j;
		int iP = P.i;
		int jP = P.j;

		/*std::cout << "(" << BLACS::myrank << ") multiplying " << A.x << "x"<<A.y<<" with a "<<B.x<<"x"<<B.y<<", alpha = "<<alpha << " opA="<<opA<<" opB="<<opB<<" blocks are "<<A.rblock()<<"x"<<A.cblock()<< " and "<<B.rblock()<<"x"<<B.cblock()<< std::endl;*/

		/*PBLAS::pxgemm(opA, opB, m, n, k, (Scalar)1, 
			A.localData(), iA, jA, A.descriptor(),		// matrix A
			B.localData(), iB, jB, B.descriptor(),	// B
			(Scalar)0, P.localData(), iP, jP, P.descriptor());	// resulting temporary product	*/
		
		//if (P.y == 1)
		//{
		//	// Unlike *gemm, for Matrix-Vector multiplication m,n must be given before any transpose operation...
		//	// x and y are swapped upon transpose(), so we have to un-swap them here
		//	PBLAS::pxgemv(opA, (A.use_transpose ? A.y : A.x), (A.use_transpose ? A.x : A.y), alpha, 
		//		A.localData(), iA, jA, A.descriptor(),
		//		B.localData(), iB, jB, B.descriptor(), /*incx*/1, 
		//		/*beta*/(Scalar)0, P.localData(), iP, jP, P.descriptor(), /*incy*/1);
		//}
		//else
		//{		

		

		/*if(A.local_matrix.size()>0)
			std::cout << "(" << BLACS::myrank << ") A(1,1)=" << *(A.localData()) << std::endl;*/

		PBLAS::pxgemm(opA, opB, A.x, B.y, A.y, alpha, 
			A.localData(), iA, jA, A.descriptor(),		// matrix A
			B.localData(), iB, jB, B.descriptor(),	// B
			(Scalar)0, P.localData(), iP, jP, P.descriptor());	// resulting temporary product
		//}
		//std::cout << "(" << BLACS::myrank << ") WTH"<< std::endl;

		A.clear();
		B.clear();
		return P;
	}


}	// end namespace peigen
#endif // PEIGEN_SHARED_PROD_H